# PLAN — Rust rewrite of `segment`

Reimplementation of Harward & Woodcock's nested-hierarchical region-growing
segmenter (`~/mac2025/segment`, C) in Rust.

**Definition of done:** byte-exact reproduction of both golden outputs — `rmap.51`
and `armap.58` for Case 1, `rmap.17` and `armap.1` for Case 2 — at latency
comparable to the C program, and able to segment a 15000 × 15000 × 6 image that
the C program cannot.

Source read in full before writing this: `src/{main,segment,region,pixel,set,gdal_io}.c`,
`src/segment.h`, `libipw/util/{allocnd,ecalloc}.c`, `libipw/pixio/pvwrite.c`,
`inc/{typedef.h,PORT/linux/config.h}`. The paper was not needed — the C is
unambiguous about every decision below.

(The paper itself is not redistributed here; it is copyrighted. See
<https://doi.org/10.1080/01431169208904109>.)

---

## 1. What the program actually does

Two phases, both in one run. The C driver is `segment()` → `main_loop()` → `wind_up()`.

**Phase 0 — pixel pre-merge** (`pixel.c: pixel_pass`)
Before any region exists, two sweeps over the raw image:
- `pix_nnbr`: for every pixel, compute squared distance to each of its 4 (or 8)
  neighbours in integer arithmetic. If the minimum is ≤ Tg², set a bit in the
  contiguity byte for *every* direction achieving that minimum.
- `pix_merge`: pair up mutually-nearest pixels, one merge per pixel, giving each
  pair (or singleton) a region id. This is what turns 62500 pixels into 55226
  regions in Case 1 — an 12% reduction whose only purpose is to shrink the
  initial region list.

**Phase 1 — normal passes** (`region.c: seg_pass`, driven by `main_loop`)
Repeat until a pass produces zero merges:
1. For every active region, find its nearest neighbour by centroid distance
   (`reg_nnbr`), recording id and d².
2. If a merge coefficient < 1 is set, histogram the d² values and pick a *pass*
   tolerance Tp² admitting at most `cm × nreg` merges (`get_tp2`).
3. For every active region in id order, merge with its nearest neighbour if:
   d² ≤ Tp², the neighbour is still active and unmerged this pass, the pairing is
   mutual (within FLT_EPSILON), at least one side is below `nviable`, and the
   union is ≤ `nmax`. Merges always fold into the **lower** id.

On termination write `<base>.rmap.<pass>`. Loop over remaining `-t` tolerances.

**Phase 2 — auxiliary passes** (`region.c: seg_apass`, driven by `wind_up`)
Repeat until no merges: force every region below `nnormin` pixels (or `nabsmin`,
for "special" regions under `-B/-N`) to merge with its nearest neighbour —
**with no distance ceiling**. Size cap relaxes from `nmax` to `nabsmax`.
On termination write `<base>.armap.<apass>`.

The bounding-box discipline is the load-bearing structure of the whole program:
no region ever stores a pixel list. `reg_nnbr`, `merge_regions`, and
`mark_reg_in_image` all scan a region's bounding box and test `regid == r`. The
bbox is the only spatial index.

**Two exits to comment out, per your instruction:**
- `segment.c:wind_up` opens with `if (Spr->nnormin == 1) exit(0);`
- `-S breakpoint` dumps `nnbrlist` to `spr_<nreg>` and exits before Phase 2.

The rewrite always runs both phases and always writes both maps. `-S` is dropped
entirely — it exists only to checkpoint between the phases.

---

## 2. Ground truth

| | Case 1 | Case 2 |
|---|---|---|
| Input | `tests/golden/misc/temp_byte_bip` (ENVI) ≡ `test_3456/input/test_3456.bip.ipw` | `LC80220492014083LGN00/input/…_stack.ipw` |
| Shape | 250 × 250 × 4, uint8 | 250 × 250 × 8, uint8 |
| Phase 1 | 51 passes → `rmap.51` | 17 passes → `rmap.17` |
| Phase 2 | 58 passes → `armap.58` | 1 pass → `armap.1` |
| Initial regions | 55226 / 62500 (0.88) | 31609 / 62500 (0.51) |

```
segment -t 10 -m .1 -n 15,15,100,2500,2500 -o t10-m1-n15_15_100_2500_2500_myseg <input>
```

⇒ `nabsmin=15, nnormin=15, nviable=100, nmax=2500, nabsmax=2500`, 4-way, no mask,
no log band, no normality band.

Comparison target is the raw payload: 125000 bytes = 250·250·2, uint16 LE.
`tests/golden/test_3456/expected/proof/regmap.armap.58` and `…rmap.51` are exactly
that. See `tests/GOLDEN.md` for the IPW-vs-ENVI container note and the
`temp_byte_bip` / `test_3456.bip` trap.

---

## 3. Bit-exactness hazards, ranked

This is where the rewrite will actually fail. Ordered by how likely each is to
cost a day.

### 3.1 `flip()` — the global RNG (highest risk)

```c
#define flip()  (random() & 01)
```

`random()` is **never seeded**, so it runs from the default seed of 1. It is
called in `reg_nnbr` whenever a candidate neighbour's distance exactly equals the
running minimum, and (only under `-A`) in `seg_apass`. Early passes over uniform
imagery produce many exact ties, so this is consumed constantly and the output
depends on the whole call sequence.

**Resolved (Q1): the golden files came off Linux, so the target is glibc.**
glibc `random()` is TYPE_3: 31-word additive feedback, `r[i] = r[i-3] + r[i-31]`,
output `(u32)r[i] >> 1`; state seeded by the Lehmer recurrence
`16807·r[i-1] mod 2^31-1` in Schrage form, then 310 outputs discarded. ~30 lines.

**Correction — macOS `random()` is NOT a different generator.** This plan
originally asserted it was, and that was wrong. Verified directly: an unseeded
`random()` on this machine emits `1804289383, 846930886, 1681692777, …` —
byte-identical to glibc. Both descend from 4.3BSD TYPE_3 and Apple's libc uses
the same 16807-Schrage seeding. The C reference reproduces the golden files with
the shim *removed*.

The port is still the right call — it pins the generator against libc version
drift and against ever building on a platform that does differ — but it was never
what stood between us and the golden bytes.

**What actually blocked the C reference on macOS was `set.c`'s undefined
behaviour** (PLAN.md section 4, bug list): `add_to_set` reads a 4-byte
`REGION_ID` through a `long *` and dedups on all 8 bytes. On Linux/x86-64 the
adjacent stack garbage happened to be stable within a `reg_nnbr` call, so it
behaved as a 32-bit compare. On macOS/arm64 it is not stable, duplicate ids
escape dedup, and each duplicate creates an exact-distance tie that consumes an
extra `flip()` draw. That is the failure this section was worried about, arriving
by a different route: not a different RNG, but a different number of draws.

**`getpagesize()` is a second host-dependent hazard.** `main.c` derives
`reclaim_trigger` from it, which sets *when* `compact_region_list` renumbers:
4096 gives 911 on Linux, 16384 gives 3641 on Apple silicon. Measured both ways —
the region-map payloads stay byte-exact either way, so compaction really is
output-neutral (id renumbering preserves ascending order), but the trigger has to
be pinned to the Linux value for a line-for-line `myseg.log` match. The Rust
hardcodes 4096; independently, the compaction pattern in the golden log pins the
trigger to (869, 939], and 911 falls inside it.

Reproducing the stream also requires the call count and order to match, which
means bit-identical f32 distances (§3.2) and identical neighbour-set order (§3.3).

**The tie-break experiment, run:** ties are load-bearing. Forcing `flip()` to a
constant diverges the output both ways, and the two constants diverge from each
other:

| build | rmap.51 | armap.58 | myseg.log diff |
|---|---|---|---|
| glibc RNG | **match** | **match** | 3 cosmetic lines |
| `flip() = 0` | differ | differ | 1300 lines |
| `flip() = 1` | differ | differ | 1094 lines |

All three still converge at 51/58 passes. Two further C-side results: `-O2` and
even `-ffp-contract=fast` stay byte-exact on this input. The latter is one input,
not a proof, so Rust keeps contraction off regardless.

### 3.2 Floating-point arithmetic

Centroids are `f32` and every distance is accumulated in `f32`:

```c
Ctr1[band] = (R1->npix * Ctr1[band] + R2->npix * Ctr2[band]) / (R1->npix + R2->npix);
```

Reproducing this needs care in three places:
- `npix` is `ushort` → promoted to `int` → converted to `float` for the products;
  the divisor `(n1 + n2)` is computed in **int** and then converted. Rust must
  mirror that, not compute the divisor in f32.
- The `reg_dist2` loop accumulates band-by-band in f32 in band order. Any
  reassociation or FMA contraction changes the result. Rust does not contract or
  reassociate float ops by default, so plain `f32` code is safe — but nothing in
  the build may enable fast-math, and no `mul_add`.
- On x86-64 (SSE2) and aarch64, C's `FLT_EVAL_METHOD` is 0, so intermediates are
  genuinely f32 and Rust matches. This would not hold on 32-bit x87.

The mutual-nearest test is `fabs(a - b) <= FLT_EPSILON` — an f32 subtraction
widened to double, compared against `FLT_EPSILON` as a double. Reproduce that
widening exactly; comparing in f32 is not the same predicate.

**Centroids cannot be recomputed from integer sums.** The repeated rounded
weighted average drifts from the true mean. This forecloses the obvious memory
optimisation in §6.

### 3.3 Neighbour set iteration order

`set.c` is an insertion-ordered array with linear backwards dedup — *not* a hash
set. `reg_nnbr` iterates it in insertion order, and that order decides which
candidate first establishes the running minimum and therefore how many `flip()`
calls follow. The Rust version must use an insertion-ordered vector with linear
dedup and populate it in the same scan order: bbox rows top-to-bottom, samples
left-to-right, directions in `cd4`/`cd8` order.

Capacity is `MAX_NEIGHBORS` = 5000; overflow is a hard error in C. Keep that.

### 3.4 Integer distance in Phase 0

`pix_dist2` accumulates in `long` from `uchar` operands — exact, no float. Use
`i64`/`u32` and it is trivially reproducible. Note the asymmetry: Phase 0 compares
**pixels** in integer space, Phase 1+ compares **centroids** in f32 space.

### 3.5 Traversal order

Every loop that can affect output is in a fixed order: pixels row-major, regions
by ascending id, directions by table index. `pix_merge` additionally carries a
rotating `idir` cursor across pixels so the direction search starts somewhere
different each time — a stateful detail that is easy to miss and changes which
pixel pairs form.

---

## 4. The C memory tricks — keep or drop

| # | Trick | Why C needed it | Verdict for Rust |
|---|---|---|---|
| 1 | **Dope vectors** (`allocnd`: array of row pointers into one block) | Gives `img[l][s]` syntax without 2-D array support | **Drop.** Costs 8 bytes/row and an indirection. Flat `Vec<T>` + `y * width + x`. Purely a C ergonomics hack. |
| 2 | **Contiguity band: 1 byte/pixel, meaning changes 3× through the program** | Genuine 8× saving vs. a byte per direction | **Keep the representation, drop the overloading.** 225 MB at 15000² is worth having. But encode the three meanings (pixel-NN flags → chosen merge direction → boundary/foreign-region mask) as three distinct newtypes over `u8` so the transitions are typed instead of documented in a 30-line comment. |
| 3 | **`ushort npix`, `short` bbox coords** (region = 12 bytes) | Halves the largest array | **Keep 12 bytes, but use `u16` coords, not `i16`.** The C limit is `MAXSHORT` = 32767 because the coords are *signed* — nothing stores a negative coordinate. Unsigned lifts the ceiling to 65535 per axis for free. `npix` stays `u16`, which is already implied by the CLI validating `nabsmax ≤ 65535`. See §11 Q5. |
| 4 | **`REGION_ID` = `unsigned` (u32)** | — | **Keep.** 4.29e9 ids covers 65535² pixels. |
| 5 | **Reusing `nnbrlist[r].nbr_id` as the old→new id table during compaction** | Avoids a second n-element array | **Keep.** At 15000² a separate table is another 800 MB. It is only an aliasing trick, not UB — express it as an explicit "the neighbour list is now a translation table" state rather than a silent reuse. |
| 6 | **`SITEM` union + runtime size dispatch in `set.c`** | Poor man's generics | **Drop.** `Vec<u32>`. The union wastes 4 bytes/entry and the `case 4` branch reads 8 bytes out of a 4-byte object — real UB the source already flags as `// OFFENDER`. It happens to work because the garbage upper half is stable within a call. |
| 7 | **Bounding-box scan instead of per-region pixel lists** | The whole point — pixel lists would dwarf everything else | **Keep.** This is algorithmic, not a C workaround. It is also the dominant cost (§7). |
| 8 | **Region 0 reserved for masked pixels; dummy region at `nreg+1`** | 1-based indexing without offsets | **Keep.** Free, and it keeps ids identical to the C. |
| 9 | **`realloc` shrink in `compact_region_list`** | Returns memory mid-run | **Keep** as `truncate` + `shrink_to_fit`. |
| 10 | `register` / `REG_1..REG_6` macros | 1989 | **Drop.** No-ops for thirty years. |

**Summary:** tricks 1, 6, 10 are pure C artifacts. Trick 5 is a C artifact that
is still worth keeping for its memory profile at scale. Tricks 2, 3, 4, 8, 9 are
real space optimisations that remain correct and worthwhile in Rust. Trick 7 is
the algorithm.

### Bugs found while reading — decide deliberately, do not port blindly

- `segment.c:~120` — `for (band = 0; band > Spr->nbands; band++)` zeroing region 0's
  centroid never executes (`>` for `<`). Harmless: `ecalloc` already zeroed it.
- `main.c` — `sproc.cm` is read (`if (sproc.cm <= 0. …)`) before being
  initialised when `-m` is absent. Reading an uninitialised stack `float`. The
  tests always pass `-m`, so it is latent. Rust: default `cm = 1.0`.
- `set.c:get_from_set` case 4 writes through `int*` what `add_to_set` read through
  `long*`. Both vanish with `Vec<u32>`.
- `merge_regions` errors out if `npix1 + npix2 > 65535`. Retain as a checked error,
  not a wrap.

---

## 5. Why the C segfaults at 15000 × 15000

Not the `short` coordinates — 15000 < 32767, so `do_headers`' `MAXSHORT` guard
passes. Three 32-bit overflows downstream:

1. `ecalloc(int nelem, int elsize)` and `allocnd`'s `int bsize` are **signed 32-bit**.
   The centroid list is `ecalloc(nreg + 2, nbands * sizeof(float))`. At 15000²
   with ~0.88 initial regions/pixel and 6 bands that is 198e6 × 24 ≈ **4.75e9**,
   which overflows `int` to a negative value → `assert(nelem > 0)` or a wild
   allocation.
2. `read_image`/`GDAL_read_image` compute `image_size = nlines * nbands * nsamps`
   in `int`: 1.35e9 for 6 bands is within range, but 8 bands at 15000² is 1.8e9
   and 20000² overflows outright.
3. `Spr->nreg` is `long` but is printed with `%d` and passed through `(int)` casts
   at both allocation sites.

Rust fixes all three by using `usize` throughout. The remaining constraint is
real memory, §6.

---

## 6. Memory budget at 15000 × 15000 × 6

225e6 pixels. Initial region count scales with content; Case 1's 0.88
regions/pixel is the pessimistic end, giving ~198e6 initial regions.

**Measured, 15000 x 15000 x 6:** analytic peak 6.45 GB, actual peak RSS
5.48 GB, wall clock 157.9 s. The table below was the pre-implementation
estimate; the initial region ratio turned out to be 0.50, not 0.88.

| Array | Element | Bytes |
|---|---|---|
| image (u8, BIP) | 6 | 1.35 GB |
| contiguity band | 1 | 0.23 GB |
| region band | 4 | 0.90 GB |
| region list | 12 | 2.38 GB |
| **centroid list** | **24** | **4.75 GB** |
| neighbour list | 8 | 1.58 GB |

The image is freed before the neighbour list is allocated, so peak ≈ **9.9 GB**,
hit at the moment the initial region list is built. It falls fast: the first
compaction cuts the region-side arrays roughly in half, and by the end of Phase 1
they are negligible.

The centroid list dominates and **cannot be shrunk without breaking exactness**
(§3.2) — storing `u32` band sums and dividing on demand would give the true mean,
which is not what the C computes.

Levers, if 9.9 GB is too much (§11 Q3):
- Build the region list in id-order chunks, spilling centroids that no active
  region references — complex, and the access pattern is not chunk-local.
- Memory-map the region and centroid lists. Keeps exactness, trades RAM for page
  faults; the bbox scan pattern has poor locality so this could be very slow.
- A `--max-regions` guard that fails fast with a clear message rather than an OOM kill.

**Resolved (Q3):** 9.9 GB is acceptable, so ship the straightforward version and
treat spilling as a separate project if it is ever needed.

The typical workload — a 5000 × 5000 × 6 tile — is 11× smaller: ~0.9 GB peak,
comfortable anywhere. Nodata reduces it further, since masked pixels never enter
the region list at all.

---

## 7. Performance

The C is single-threaded, so parity is the easy part; 15000² is where it matters.

**Dominant cost:** `reg_nnbr` runs once per active region per pass and scans the
region's full bounding box, testing `regid == r` on every pixel. A long thin
region wastes most of that scan. With 51 + 58 passes, this is essentially all the
runtime. `merge_regions` scans the union bbox twice more.

**Plan:**
1. Get it exact and single-threaded first. Measure against the C on the 250²
   cases and on a synthetic 5000².
2. Cheap wins that cannot change results: flat indexing over dope vectors, row
   slices instead of per-pixel bounds checks, `SmallVec` for the neighbour set to
   avoid touching the 5000-entry buffer, skipping `Cinternal` pixels early
   (already in the C).
3. **Parallelising `reg_nnbr` without breaking the RNG.** The first loop of
   `seg_pass` is read-only over `rband`/`cband`/`ctrlist`, so it parallelises —
   except that `flip()` is a serial global. Split it:
   - *Parallel:* per region, scan the bbox and produce the neighbour candidates
     in insertion order with their f32 distances.
   - *Serial, in ascending id:* replay the exact C selection loop over that list.
     A `flip()` is consumed precisely when a candidate's distance equals the
     running minimum, so replaying the ordered list reproduces the RNG stream
     call-for-call.

   Process in id-ordered chunks to bound the memory held by the parallel stage.
   This keeps byte-exactness while moving the expensive part off the critical path.
4. The merge loop stays serial — it is order-dependent by construction.

---

## 8. Rust structure

```
src/
  main.rs        CLI (clap), argument validation mirroring main.c's rules
  config.rs      SegConfig: tolerances, cm, the five -n values, flags
  image.rs       Image { data: Vec<u8>, nlines, nsamps, nbands }, BIP, flat indexed
  contig.rs      The three contiguity-byte phases as distinct newtypes
  pixel.rs       pix_nnbr, pix_merge, make_region_list, pix_check_bounds_and_mask
  region.rs      RegionList (SoA), reg_nnbr, merge_regions, compact, seg_pass, seg_apass
  nbrset.rs      Insertion-ordered dedup set, capacity 5000
  mask.rs        Explicit -M mask + derived nodata mask, combined (§9.1)
  rng.rs         glibc TYPE_3 random() port; never calls platform random()
  io/
    mod.rs       Format sniffing by extension + magic
    envi.rs      ENVI raw + .hdr sidecar (read and write)
    ipw.rs       IPW header parse/emit, byte-aligned pixels, nbits masking
    tiff.rs      via `tiff` crate
    png.rs       via `png` crate
tests/
  golden.rs      Runs both cases end to end, compares against tests/golden/
```

**Region storage:** structure-of-arrays, not the C's array-of-structs —
`Vec<BBox>`, `Vec<u16>` npix, `Vec<u8>` flags, `Vec<f32>` centroids (nbands-strided).
The hot loops touch bbox and centroid but rarely flags, so SoA improves locality
without changing a single result.

**Crates:** `clap`, `tiff`, `png`, `rayon` (step 7.3 only), `smallvec`. No `gdal`
by default — the original's GDAL dependency covers exactly two formats we can
write by hand, and linking libgdal is a large cost for that. See §11 Q2.

---

## 9. I/O, mask, and nodata

**Read:** ENVI (raw + `.hdr`), IPW, TIFF/GeoTIFF, PNG. All present internally as
uint8 BIP. The original rejects anything but Byte and so does the rewrite —
Case 2's ENVI `_stack` is int16 and is *supposed* to fail; only its `.ipw` works.

IPW headers are plain-text records terminated by a form feed; pixels are
byte-aligned at `bytes` per pixel and masked to `bits` — **not** bit-packed.
`nbands` comes from the `basic_image_i` record, per-band width from `basic_image N`.

**Write:** the region map, sized to the smallest type holding `nreg`
(`nbits ≤ 8` → u8, `≤ 16` → u16, else u32), exactly as `GDAL_write_image` does.
Default output is ENVI raw + `.hdr` so the bytes line up with `proof/`; IPW output
is available for regenerating `.armap.58`-style containers. Geotransform and
projection pass through where the format carries them.

### 9.1 Nodata (water, non-treed area)

The C already has the machinery: a masked pixel gets `REGION_ID 0`, is skipped by
`pix_nnbr` and `pix_merge`, is skipped by `make_region_list`, and
`pix_check_bounds_and_mask` sets its neighbours' contiguity bits so nothing ever
tries to merge across it. Region 0's centroid stays at zero and is never used.

So nodata does not need new algorithm — it needs to be **funnelled into the
existing mask**:

```
effective_mask[p] = 0  if  p is nodata            (derived, see below)
                 or  explicit -M mask says 0
                 else 1
```

Nodata is derived from, in precedence order:
1. `--nodata <value>` on the command line;
2. the format's own declaration — ENVI `data ignore value`, GeoTIFF
   `GDAL_NODATA`, PNG `tRNS`;
3. nothing, if neither is present.

**Multi-band rule.** A pixel counts as nodata when **all** bands equal the nodata
value (`--nodata-any` switches to "any band"). All-bands is the safe default for
this use: masked-to-land imagery carries 0 across every band over water, while a
legitimate 0 in a single band is ordinary dark ground that should still segment.

**Output.** Masked pixels are written as region 0 and the output header declares
`data ignore value = 0`, so water/non-treed area round-trips as nodata rather than
as a spurious stand.

**Correctness risks specific to this path**, to be covered by tests:
- A nodata pixel must never contribute to any centroid. It has no region, so it
  cannot — but `pix_nnbr` must skip it *before* computing distances, or masked
  values leak into the Phase 0 tie structure.
- A region adjacent to nodata must not treat the nodata side as a neighbour with
  distance 0. `pix_check_bounds_and_mask` prevents this by marking the direction
  contiguous; getting this wrong produces regions that grow along shorelines.
- An image that is *entirely* nodata, and a nodata region that splits the image
  into disconnected components, should both terminate cleanly.

---

## 10. Milestones

Checked off as they land. Each gate is a command whose output is pasted into the
commit message.

- [x] **M0 — C reference on this machine.** Build the original; add
      `glibc_random.c` (§3.1); regenerate Case 1 from `misc/temp_byte_bip`.
      *Gate met.* The C in `reference/csegment/` reproduces both golden files
      byte-exactly on macOS/arm64 in 0.56 s, with `myseg.log` matching on every
      numeric value across all 51 + 58 passes (3 cosmetic tab/space lines differ).
      The blocker was `set.c` UB, not the RNG — see section 3.1. Tie-break
      experiment run; ties matter.
- [x] **M1a — ENVI + IPW.** Readers and the ENVI region-map writer.
      *Gate met:* `tests/io_golden.rs`, 5 tests. ENVI and IPW readers agree
      byte-for-byte on Case 1; Case 2 IPW reads as 250x250x8; int16 ENVI is
      rejected as the original does; the writer reproduces
      `proof/regmap.armap.58` byte-exactly from its own region ids.
- [ ] **M1b — TIFF + PNG readers, IPW writer.** Not on the critical path;
      neither test case needs them.
- [x] **M2 — Phase 0.** `pix_nnbr`, `pix_merge`, `make_region_list`.
      *Gate:* `nreg` = 55226 (Case 1), 31609 (Case 2) — the numbers `myseg.log`
      records at "of a possible 62500 regions are required".
- [x] **M3 — Phase 1.** `reg_nnbr`, `merge_regions`, `seg_pass`, `compact`.
      *Gate met.* Every algorithmic line of both `myseg.log` files reproduced
      (450 lines Case 1, all counters, all passes); the only diff is the C's
      predicted-malloc-sizes block. `rmap.51` and `rmap.17` byte-exact.
- [x] **M4 — Phase 2.** `seg_apass` and `wind_up` without the two `exit(0)`s.
      *Gate met.* `armap.58` and `armap.1` byte-exact. **Definition of done
      reached** — and it retroactively proves the glibc `random()` port, since a
      single desynced draw would diverge the map.
- [x] **M5 — Mask and nodata.** `-M`, `--nodata`, `--nodata-any`, derived-mask
      plumbing (section 9.1). *Gate met:* `tests/nodata.rs`, 7 tests — masked
      pixels are region 0, a region never grows across nodata (with an unmasked
      control proving the test discriminates), nodata never reaches a centroid,
      derived nodata matches an explicit mask, an all-nodata scene terminates,
      and the multi-band all-bands rule behaves.

      *Note on a test that was NOT written:* "a nodata border segments like the
      cropped image" is **false**, and faithfully so. `pix_check_bounds_and_mask`
      sets three bits per image edge (`N_EDGE` = NW|N|NE) but only the one
      matching direction per nodata neighbour, so an edge pixel and a
      nodata-adjacent pixel can differ on the `== Cinternal` test. The C has the
      same asymmetry; it is not a bug to fix.
- [x] **M6 — Scale.** *Gate met.* 15000 × 15000 × 6 (1.35 GB input, 225M pixels)
      completes in **157.9 s** with **5.48 GB peak RSS**, writing both maps at
      900 MB each. 23 normal passes, 1 auxiliary. This is the size that segfaults
      the C. 5000 × 5000 × 6 — the stated typical tile — runs in **16.0 s** at
      **720 MB**.

      Measured peak came in well under the section 6 estimate of 9.9 GB: the
      initial region ratio on real-ish imagery is 0.50 regions/pixel, not the
      0.88 worst case, and the region arrays halve after the first pass, so RSS
      never touches the analytic peak. `--mem-report` output is printed at the
      start of every run.
- [x] **M7a — Single-threaded performance.** C -O0 24.49 s, C -O2 12.63 s, Rust 13.30 s on 5000^2 x 6 (mean of 3). Rust is 1.8x faster than the C as built and 5% slower than an -O2 C. At 5000^2 the Rust and C outputs are byte-identical (100 MB each), 400x the test-case area.
- [x] **M7b — Parallel `reg_nnbr`** (section 7.3). *Gate met:* byte-exact and
      roughly 2x faster.

      | scene | serial | parallel (10 cores) |
      |---|---|---|
      | 5000^2 x 6 | 13.3 s | **9.72 s** (mean of 3) |
      | 15000^2 x 6 | 157.9 s | **77.6 s**, 4.97 GB peak |

      Against the C on 5000^2 x 6: **1.3x faster than `-O2`**, 2.5x faster than
      the `-O0` build the original Makefile produces. Output at both sizes is
      byte-identical to the serial run, and at 5000^2 still byte-identical to the
      C.

      The split is the one section 7.3 proposed. The C's first loop only *reads*
      the region band, contiguity band and centroids; the sole order-dependent
      thing in it is the `flip()` stream. So the bbox scan fans out across
      threads, and selection replays serially in ascending id order — a draw is
      consumed exactly when a candidate ties the running minimum, so replaying
      the ordered candidate list reproduces the stream call for call. Chunked at
      2^18 regions because holding candidate lists for all 113M regions at once
      would be gigabytes.

      Speedup is ~2x rather than ~10x because only the scan parallelises: the
      selection replay and the whole merge loop stay serial by construction.
      Beating that means giving up byte-exactness, which the tests forbid.

      Threshold is 200k regions (below it, fan-out costs more than the scan);
      `--threads 1` forces serial, and `tests/segment_golden.rs` runs the golden
      cases through the parallel path with the threshold forced to 0.

M3's gate is the one that will actually save time. `myseg.log` carries `nreg`,
`dmin2`, `maxpix` and seven merge counters for all 51 passes; diffing that against
our own log localises a divergence to one pass and one rejection reason, instead
of leaving 125000 mismatched bytes and no idea why.

---

## 11. Decisions taken

| | Decision |
|---|---|
| **Q1 RNG** | Golden is Linux/glibc. Port glibc TYPE_3 `random()`; shim it into the C reference too, since macOS `random()` differs (§3.1). |
| **Q2 Formats** | ENVI, IPW, TIFF, PNG. Native readers/writers, no libgdal. uint8-only input, matching the original — int16 input is rejected, not converted. |
| **Q3 Memory** | 9.9 GB peak at 15000²×6 is acceptable. Typical workload is 5000² tiles (~1.1 GB peak), so the common case is comfortable. |
| **Q4 Exactness** | Byte-exact on the 250² test cases. 15000² is a latency target with no oracle. Exact-by-default stays; a `--fast` mode (parallel merge, deterministic tie-break in place of `flip()`) becomes legitimate for large images and is deferred to M7. |
| **Q5 Region cap** | Keep `u16 npix` — 65535 pixels/region, region struct at 12 bytes. |
| **Q6 Flags** | In: `-t -m -n -o -8 -M -b -l -B -N -A`, plus `--nodata`/`--nodata-any`. Out: `-h` (hsegment debug files), `-S` (phase-1 breakpoint). Both `exit(0)`s in `wind_up` removed — both phases always run, both maps always written. |

### Still open

**Nodata multi-band rule.** Defaulting to "all bands equal the nodata value"
(§9.1), with `--nodata-any` for the other reading. If your imagery marks nodata in
a single band — or uses a per-band value rather than one shared sentinel — tell me
and I will flip the default.

**Nodata and `nabsmin`.** Regions touching a nodata boundary can end up below
`nnormin` with no valid neighbour to merge into, since Phase 2 cannot merge across
nodata. The C would leave them and report them under "WARNING! Questionable
regions". I plan to keep that behaviour rather than invent a rule. Say so if you
want small shoreline stands handled differently.

---

## 12. Modernisation (post-M7)

The port is finished and byte-exact; these are deliberate departures from 1992
behaviour, taken one at a time with the golden check green after each.

- [x] **12.1 -- 16-bit input.** The original accepted uint8 only. Modern
      multispectral imagery is 12-bit in a 16-bit container, so an 8-bit-only
      reader forces a rescale that changes the segmentation before it starts.
      Input now carries `u8`, `u16` or `i16` (`Samples` in `image.rs`), int16
      being the container Landsat Collection 2 surface reflectance ships in.

      *Scope.* Nothing downstream of Phase 0 sees a pixel -- centroids are f32
      and the image is freed once the region list exists -- so widening is
      confined to `image.rs`, `pixel.rs` and `RegionList::from_pixel`. Phase 0
      dispatches once on the sample type and is monomorphic below that, so the
      8-bit path compiles to what it did before.

      *On f32.* An earlier note here said the wide path would need f64
      distances. On re-examination it does not, and the reason matters: the
      precision that counts is precision *near the tolerance*, and near the
      tolerance the magnitudes are small. `RegionList::dist2` only loses
      absolute precision at distances far larger than any threshold, where the
      comparison is already unambiguous. Centroids are exact in f32 for every
      `u16` and `i16` (65535 < 2^24). The one place the width was genuinely
      marginal is `pix_nnbr`'s `mdist2 <= tg2`, which now compares in f64 --
      and for 8-bit input that is provably not a change, since `mdist2` tops
      out at `255^2 * nbands`, exact in both widths. So there is one code path,
      not two.

      *Verified.* Golden: both cases, both phases, payloads byte-identical
      (125000/125000 bytes each) and both `myseg.log` files matching line for
      line. `tests/wide_input.rs` runs Case 1 again with every sample widened
      to `u16` and to `i16` and requires identical maps. The fixture
      `LC80220492014083LGN00_stack` turns out to be the real int16 reflectance
      (DN 0..8990) behind the 8-bit `.ipw` -- it now reads, segments, and gives
      a different map from its own rescaling, which is the point. Timing on
      3000^2 x 6: 38.9 s vs 38.7 s baseline, i.e. unchanged.

      *Consequence.* Tolerance is in DN and does not carry across widths.
      `-t 10` on the 8-bit rescaling is about `-t 350` on the 16-bit original.

- [x] **12.2 -- `u32 npix` and a growable neighbour set.** Two 1992 type
      widths that had become limits on the answer rather than on the machine.

      *`npix`.* Was `unsigned short`, and the "no limit" settings for
      `-n Nviable,Nmax,Nabsmax` were spelled 65535 for that reason -- so a
      default run silently stopped growing a stand at 65535 pixels, a 256 m
      square at 1 m resolution. `npix` is now `u32`, "no limit" now means
      unlimited (`MAX_REGION_PIXELS`), and `-n` accepts values above 65535.
      The old ceiling is still available by asking for it explicitly, and
      `tests/region_size.rs` pins both directions: a uniform 300x300 field now
      ends as one region of 90000 pixels, and `-n ...,65535,65535` still caps
      it. Cost is 2 bytes per region, 226 MB at 15000^2 against a 5 GB peak.

      *Neighbour set.* `MAX_NEIGHBORS = 5000` aborted the whole run on the
      5001st neighbour. The list now grows. Because an unbounded backwards
      linear dedup is quadratic, past `LINEAR_LIMIT = 96` entries membership
      moves to a boxed hash set -- `items` stays in the same insertion order
      either way, which is the property the RNG replay depends on. Verified
      against the previous commit: a 2501-pixel row flanked by 5002 singleton
      regions fails there with `more than 5000 neighbors of region 2502` and
      completes here.

      *Cost.* 3000^2 x 6: 38.4 s against 37.6 s before, about 2%. The first
      cut of the neighbour set cost 12% -- an inline `HashSet` and a split
      distance array, in a struct that exists once per scratch slot in the
      hot sweep. Boxing the side-table and keeping the `(id, dist)` pairs in
      one vec brought it back. Golden output unchanged throughout.
- [x] **12.3 -- Provenance in output headers.** IPW recorded
      `history = segment -t 10 -m .1 -n ... ../LC80220492014083LGN00_stack.ipw`
      in every image it wrote, and that record is the only reason the
      invocation behind the golden fixtures was recoverable eleven years later.
      Our ENVI output carried nothing.

      Every written map now carries the command that made it: ENVI as
      `history = {...}` plus `software = {...}`, TIFF as `ImageDescription`
      and `Software`. Arguments are shell-quoted, and `{`/`}`/newlines are
      replaced rather than escaped -- a mangled history line is better than an
      unparseable header, and `tests/provenance.rs` pins that a path
      containing a brace still leaves the header readable by our own parser.

      No timestamp, deliberately: the same command twice produces identical
      files, which is worth more here than knowing the hour of the run. The
      raster is untouched, so the golden payload comparison is unaffected --
      re-verified byte for byte after the change.
- [x] **12.4 -- `-b`/`-l`/`-B`/`-N`/`-A`.** Half-present is worse than either
      state, so each was taken to one end.

      *Wired up.* `-B band` and `-N low,high` (normality band and interval):
      a region whose centroid in that band falls outside the interval is
      *special* and is held to `Nabsmin` rather than `Nnormin` in Phase 2. The
      logic was already there and correct; there was no way to reach it. Both
      are required together, as in the C, and the C's `high <= 255` becomes a
      check against the input's actual sample range. The two extra auxiliary
      log lines the C prints under `-B` are printed too.

      `-A` allocated the mask, filled it in during Phase 2 and then dropped it
      on the floor. It is now written as `<base>.armask.<pass>`, one uint8
      band, the way `wr_armm` did.

      *Deleted.* `log_band` (`-b`), `lthr` and `lincr` (`-l`) drove the
      per-pass single-band `.log.<n>` debug files -- the same category as
      `-h`, which decision Q6 already dropped. Nothing read the fields. They
      are gone rather than left looking like features.

      `tests/flags.rs` pins the behaviour rather than the plumbing: a bright
      field with 25 isolated dark specks loses them to Phase 2 without
      `-B`/`-N` and keeps them with, and an interval that covers every centroid
      reproduces the no-`-B` run exactly.

---

## 13. Two-input segmentation (Ye et al. 2025)

The 1992 algorithm segments one image. The modification this repo now has to
support segments **two**: Woodcock & Harward's micro-segmentation on Landsat
spectral proxies as before, then a *segment-development* phase that merges those
micro-segments using a **different image over the same grid** — forest structure,
age, or species. Published as Ye, Coops, Wulder & Hermosilla, *ISPRS J.
Photogramm. Remote Sens.* 226 (2025) 381–395; the local copy is in `paper/`,
which is gitignored.

The requirement is that this stays one program. With no second image the
behaviour is exactly what it is today, byte for byte, golden fixtures and all.
Supply a second image and the auxiliary phase is replaced by the
segment-development phase.

### 13.1 What the second phase actually is

Not a re-run of Phase 2 on different pixels — the rules differ:

| | Phase 2 (1992, `seg_apass`) | Segment development (2025) |
|---|---|---|
| Input | the same image as Phase 1 | a second image, different bands, same grid |
| Starting map | the `armap` in progress | the **`rmap`** — normal passes only |
| Merge partner | mutual nearest neighbour | undersized region's nearest neighbour, **made** mutual by write-back |
| Surviving id | the lower id | the **absorbing (smaller)** region's id |
| Distance | f32 | f64 |
| Tie-break | `flip()`, glibc `random()` | none needed (see 13.3) |
| Size rules | `Nabsmin`/`Nnormin`/`Nviable`/`Nmax`/`Nabsmax`, normality band | `Nmin` and `Nmax` only |
| Masking | from the input image's nodata | **from the second image**: a region more than half non-treed is dropped |

The write-back is the substantive idea. §4.2 of the paper: *"for any segment A
that was smaller than the minimum region size, it could be merged with its
nearest neighbour segment B, even if segment A was not the nearest neighbour of
segment B."* Implemented by having A stamp its own distance onto B, so the
mutual-nearest test that follows passes by construction unless a closer A′
claims B first.

### 13.2 The oracle, and the two bugs in it

Elaine's Python (`~/mac2025/segment_python`, commit `427a5a3`) is the definition
of the phase — it produced the published results. It is vendored into
`tools/stage2_oracle/` with a reviewable diff, and `tests/stage2/` holds six
cases generated from it. `tests/STAGE2.md` is the full description; the
generator refuses to write a case it cannot pin.

Two defects had to be dealt with rather than ported:

**The whole-map wipe.** `region_map[tuple(np.array([...]).T.tolist())] = 0` with
an empty list is `region_map[()] = 0` — numpy's spelling of *the entire array*.
Whenever no region was majority-nodata, the phase zeroed every id, found no
adjacencies, merged nothing, and wrote the input back out. It never fired on the
real runs because the structure/age/species layers always contain some
majority-nodata region, so the published results are unaffected; it fires
immediately on clean input, which is why Elaine's own `integ_test.py` was
returning its 5×5 input unchanged and looking like it passed. Guarded in the
vendored copy. **Not** to be reproduced in Rust.

**The unseeded tie-break.** Near-equal candidates are chosen between with
`randint(0, 1)` while iterating a Python `set`, so the phase has no defined
answer. Sweeping all six combinations of (set/ascending/descending order) ×
(keep/take on tie), three of the six fixture cases come out different. Pretending
otherwise would have frozen one arbitrary sample as the oracle. The fixtures pin
a rule instead — **ascending region id, keep the incumbent on a near-tie** —
under which the coin flip is unreachable, so stage 2 needs no RNG at all.

Be exact about the justification, because the 1992 C has two mechanisms that are
easy to conflate. Its *merge survivor* rule is deterministic and id-based:
`if (r < nnbr_id) merge_regions(Spr, r, nnbr_id); else merge_regions(Spr,
nnbr_id, r)`, and `region.c` says so — *"merge the two regions into the region
with the lower REGION_ID"*. Its *candidate tie-break*, in `reg_nnbr`, is
`flip()` — the same unseeded coin the Python uses, commented *"This is biased,
but it does give some randomness to nnbr selection"*. So the chosen rule is
**consistent in spirit with the C's survivor convention, not a reconstruction of
its tie-break**; the C has no deterministic tie-break to reconstruct. The three
tie-insensitive cases would pass under any rule; the three sensitive ones are
what test that this one is implemented.

The oracle enforces it structurally rather than by patching from outside:
`region.py` exposes `ON_TIE` (`keep`/`take`/`random`) and only `random` reaches
`randint`. Regenerating all six cases under the default produced bytes identical
to the committed manifest, which is the proof the flip was never consulted.

A third quirk is faithful and stays: a region's centroid is the mean over **all**
its pixels including nodata ones, so non-treed pixels drag the mean toward zero.
Only the >50 % test excludes them. That is what the published segmentation did.

### 13.3 Specification

Per pass, until a pass merges nothing (region ids are visited in
**raster-first-occurrence order** throughout — the order ids first appear
scanning the region map row by row):

1. Reset every region's recorded nearest distance to +∞.
2. For each region with `npix < Nmin`: scan its bounding box; for each of its
   pixels take the 4-neighbours that are in bounds, unmasked, and a different
   region; among those distinct ids pick the smallest squared Euclidean distance
   between stage-2 centroids, **ascending by id, keeping the incumbent** when
   `|d − dbest| ≤ 1e-6·max(|d|,|dbest|)`. Record `(id, d)`. If `d` is less than
   that neighbour's own recorded distance, write `(this region, d)` onto the
   neighbour.
3. For each region with `npix < Nmin`, in the same order: skip if it or its
   partner already merged this pass, if the partner is 0, if the distance is
   still ∞, if `npix(a) + npix(b) > Nmax`, or if the two recorded distances
   differ by more than `1e-9` relative. Otherwise **a absorbs b, keeping a's id**;
   centroid becomes `(na·ca + nb·cb) / (na + nb)` in f64, in that association.
4. Rewrite the region map.

Before the first pass, centroids are per-band means of the second image over each
region, and every region whose count of all-bands-zero pixels is **strictly more
than half** its area is deleted and its pixels set to region 0.

**On arithmetic.** f64, unlike stage 1's f32 — this is not the C being ported.
The initial means are exactly reproducible without imitating numpy: the samples
are integers, every partial sum is an exact integer below 2^53, so summation
order is irrelevant and one `i64` sum plus one `f64` divide matches bit for bit.

**On adjacency.** The Python maintains a 4-bit-per-pixel adjacency band
incrementally, which forces it to keep a coordinate list per region and makes
merging quadratic. It is equivalent to recomputing the band from the region map
at the start of each pass, because regions participate in at most one merge per
pass so the updated pairs are disjoint. *Verified*: on all six fixtures the
recomputed band is bit-identical to the incremental one at every pass and the
final maps agree. So the Rust keeps no coordinate lists — bbox, `npix`, centroid
and the existing `rband` are enough, and neighbour collection reuses the shape of
`collect_nbrs`/`select_nnbr` already in `segment.rs`.

### 13.4 Shape in this codebase

- `src/stage2.rs`, new. Does **not** extend `Segmenter`: different centroid
  width, different merge direction, no RNG, no contiguity band. Sharing it would
  put a mode flag through the hot loop of a phase that is currently byte-exact.
- `io::envi::read_region_map`, new and separate from `io::envi::read`. Accepts
  ENVI data types 1, 12 and 13 — the fixtures need uint16 (our own writer) and
  uint32 (the NTEMS run). The image reader keeps refusing anything but 8- and
  16-bit; a region map is labels, not DN, and the two should not share a path.
- CLI, additive:

      --stage2 <IMAGE>     second-stage image; enables the phase
      --n2 <NMIN,NMAX>     minimum and maximum region size for it
      --rmap <FILE>        start from this region map, skip stage 1

  Absent `--stage2`, nothing changes. `--rmap` is what lets stage 2 run against
  a region map this program did not produce, which is how four of the six
  fixtures work; it also makes the phase usable on the published NTEMS stage-1
  outputs directly.
- Output stays `<base>.armap.<passes>`, with `passes` counted the Python's way
  (merging passes plus the final no-op one).

### 13.5 Milestones

- [x] **M8a — Test data.** Six cases in `tests/stage2/`, checksum-pinned and
      locked (`tests/verify_stage2.sh`, 43 files), generated from the vendored
      oracle with the arbitrariness sweep as a gate. Coverage across the set:
      every rejection branch (`over_max` 10 715 in `age_capped`, `not_mutual`,
      `busy`, `no_cand`, `inf`), single- and six-band second images, uint16 and
      uint32 region maps, a pre-masked stage-1 map, masking introduced by
      stage 2, and two cases whose stage-1 map this repo's own binary produced.
      *Gate met:* `tests/stage2_fixtures.rs`, 5 tests, deriving the invariants
      from the bytes rather than from the Python — no region split, no invented
      id, masking only grows, nothing over `Nmax` that did not start there.
- [x] **M8b — `read_region_map` + `--rmap`.** `io::envi::read_region_map` reads
      ENVI data types 1, 12 and 13 into a `u32` band, refusing multi-band and
      big-endian maps. *Gate met:* every fixture's region map is read at its
      declared type and count (`tests/stage2_fixtures.rs`), and `--rmap` drives
      a whole fixture through the command line (`tests/stage2_cli.rs`).
- [x] **M8c — The phase.** `src/stage2.rs`. *Gate met:* all six
      `expected/armap.<n>` reproduced **byte-exactly** at the oracle's own pass
      counts — 100, 250 000, 250 000, 250 000, 80 000, 80 000 bytes — plus every
      per-pass counter (`considered`, `no_cand`, `busy`, `inf`, `over_max`,
      `not_mutual`, `merged`) matching `case.json` across all ~300 passes.
      `tests/stage2_match.rs`, which also carries a negative control: perturbing
      `Nmin` or `Nmax` must move the bytes, so the gate cannot pass by being
      inert.
- [x] **M8d — Composition.** `--stage2` is threaded through
      `driver::run_with_stage2`; `run` delegates to it with `None`, so the
      one-image path is not merely equivalent to what it was, it is the same
      code. *Gate met:* both golden cases still reproduce byte-exactly
      (`tests/segment_golden.rs` unchanged), and one invocation equals two — a
      composed `-t … --stage2 …` run produces the same stage-1 map as the plain
      run and the same final map as running `--rmap` on that map
      (`tests/stage2_cli.rs`).

      The gate as originally written — "reproduces `e2e_gsv` from the proxies
      crop in a single invocation" — is also met, run by hand against the
      external drive:

          fast_segment -t 50 -m 0.2 -n 9,18,36 \
              --stage2 gsv_crop --n2 50,8000 -o tile399 proxies_crop

      writes `tile399.rmap.41`, byte-identical to `e2e_gsv/input/rmap`, and
      `tile399.armap.39`, byte-identical to `e2e_gsv/expected/armap.39`
      (80 000 bytes each). It is not in `cargo test` because the proxies crop
      lives on the drive and is not vendored — only the stage-1 map it produced
      is. The composition equivalence above is the CI stand-in.

### 13.6 Measured on the full tile

Tile 399 at 5000 × 5000 × 6 (proxies) + 5000 × 5000 × 1 (`elev_p95`),
`-t 50 -m 0.2 -n 9,18,36 --n2 80,8000`, one command: **35.6 s, 1.20 GB peak**.
Stage 1 alone with its own auxiliary phase is 31.0 s / 1.17 GB, and segment
development alone from a saved `.rmap` is 11.5 s / 0.53 GB. 4 377 977
micro-segments to 122 552 stands in 114 passes, 1 863 704 regions dropped as
majority non-treed.

No oracle exists at that size, so correctness there rests on two things:

- **A fresh 1000 × 1000 case cross-checked against the Python** — 16× the area of
  any fixture, 167 293 input regions, not part of the pinned set. Byte-identical,
  77 passes both. Rust 0.28 s, Python 14.5 s.
- **Invariants re-derived from the 100 MB output** with numpy, independent of the
  implementation: 0 pixels whose stage-1 region was split, 0 invented output ids,
  0 masked pixels that came back, 0 regions that grew past `Nmax`.

30.7 % of output regions are still under `Nmin` at 5000², and 33.0 % are in the
*Python's own* 1000² output — so that is the algorithm converging (a region with
no unmerged neighbour is stuck), not a port defect.

The `find_nearest`/`relabel` bounding-box scans are the obvious thing to attack
if this ever needs to be faster; they are O(bbox), not O(pixels in region), and a
long thin region pays for its whole box. It has not been worth it: the phase is a
third of a run that is dominated by stage 1.

### 13.7 Two counter-fidelity hazards

Both produce the correct map and the wrong per-pass numbers, which is the worst
kind of divergence: the answer looks right and the debugging aid lies.

1. **A region absorbed earlier in the same pass is still visited.** The oracle
   deletes absorbed regions from its dict only at end of pass, so the merge loop
   still reaches them, counts them in `considered`, and rejects them at the
   `busy` test. Skipping them on an `alive` flag gives the same merges and
   understates both counters.

2. **`nearest_region_id` is not reset between passes — only the distance is.**
   `update_nearest_region_dist(inf)` resets the distance; the id stands from the
   previous pass. So a region that finds no candidate this pass still carries a
   stale partner id, misses the `no_cand` branch, and lands in `inf`. Following
   a stale id is always safe: a distance is finite only if it was written this
   pass, and both writers set the id along with it.

### Still open

**Nodata multi-band rule.** Defaulting to "all bands equal the nodata value"
(§9.1), with `--nodata-any` for the other reading. If your imagery marks nodata in
a single band — or uses a per-band value rather than one shared sentinel — tell me
and I will flip the default.

**Nodata and `nabsmin`.** Regions touching a nodata boundary can end up below
`nnormin` with no valid neighbour to merge into, since Phase 2 cannot merge across
nodata. The C would leave them and report them under "WARNING! Questionable
regions". I plan to keep that behaviour rather than invent a rule. Say so if you
want small shoreline stands handled differently.

---

## 12. Modernisation (post-M7)

The port is finished and byte-exact; these are deliberate departures from 1992
behaviour, taken one at a time with the golden check green after each.

- [x] **12.1 -- 16-bit input.** The original accepted uint8 only. Modern
      multispectral imagery is 12-bit in a 16-bit container, so an 8-bit-only
      reader forces a rescale that changes the segmentation before it starts.
      Input now carries `u8`, `u16` or `i16` (`Samples` in `image.rs`), int16
      being the container Landsat Collection 2 surface reflectance ships in.

      *Scope.* Nothing downstream of Phase 0 sees a pixel -- centroids are f32
      and the image is freed once the region list exists -- so widening is
      confined to `image.rs`, `pixel.rs` and `RegionList::from_pixel`. Phase 0
      dispatches once on the sample type and is monomorphic below that, so the
      8-bit path compiles to what it did before.

      *On f32.* An earlier note here said the wide path would need f64
      distances. On re-examination it does not, and the reason matters: the
      precision that counts is precision *near the tolerance*, and near the
      tolerance the magnitudes are small. `RegionList::dist2` only loses
      absolute precision at distances far larger than any threshold, where the
      comparison is already unambiguous. Centroids are exact in f32 for every
      `u16` and `i16` (65535 < 2^24). The one place the width was genuinely
      marginal is `pix_nnbr`'s `mdist2 <= tg2`, which now compares in f64 --
      and for 8-bit input that is provably not a change, since `mdist2` tops
      out at `255^2 * nbands`, exact in both widths. So there is one code path,
      not two.

      *Verified.* Golden: both cases, both phases, payloads byte-identical
      (125000/125000 bytes each) and both `myseg.log` files matching line for
      line. `tests/wide_input.rs` runs Case 1 again with every sample widened
      to `u16` and to `i16` and requires identical maps. The fixture
      `LC80220492014083LGN00_stack` turns out to be the real int16 reflectance
      (DN 0..8990) behind the 8-bit `.ipw` -- it now reads, segments, and gives
      a different map from its own rescaling, which is the point. Timing on
      3000^2 x 6: 38.9 s vs 38.7 s baseline, i.e. unchanged.

      *Consequence.* Tolerance is in DN and does not carry across widths.
      `-t 10` on the 8-bit rescaling is about `-t 350` on the 16-bit original.

- [x] **12.2 -- `u32 npix` and a growable neighbour set.** Two 1992 type
      widths that had become limits on the answer rather than on the machine.

      *`npix`.* Was `unsigned short`, and the "no limit" settings for
      `-n Nviable,Nmax,Nabsmax` were spelled 65535 for that reason -- so a
      default run silently stopped growing a stand at 65535 pixels, a 256 m
      square at 1 m resolution. `npix` is now `u32`, "no limit" now means
      unlimited (`MAX_REGION_PIXELS`), and `-n` accepts values above 65535.
      The old ceiling is still available by asking for it explicitly, and
      `tests/region_size.rs` pins both directions: a uniform 300x300 field now
      ends as one region of 90000 pixels, and `-n ...,65535,65535` still caps
      it. Cost is 2 bytes per region, 226 MB at 15000^2 against a 5 GB peak.

      *Neighbour set.* `MAX_NEIGHBORS = 5000` aborted the whole run on the
      5001st neighbour. The list now grows. Because an unbounded backwards
      linear dedup is quadratic, past `LINEAR_LIMIT = 96` entries membership
      moves to a boxed hash set -- `items` stays in the same insertion order
      either way, which is the property the RNG replay depends on. Verified
      against the previous commit: a 2501-pixel row flanked by 5002 singleton
      regions fails there with `more than 5000 neighbors of region 2502` and
      completes here.

      *Cost.* 3000^2 x 6: 38.4 s against 37.6 s before, about 2%. The first
      cut of the neighbour set cost 12% -- an inline `HashSet` and a split
      distance array, in a struct that exists once per scratch slot in the
      hot sweep. Boxing the side-table and keeping the `(id, dist)` pairs in
      one vec brought it back. Golden output unchanged throughout.
- [x] **12.3 -- Provenance in output headers.** IPW recorded
      `history = segment -t 10 -m .1 -n ... ../LC80220492014083LGN00_stack.ipw`
      in every image it wrote, and that record is the only reason the
      invocation behind the golden fixtures was recoverable eleven years later.
      Our ENVI output carried nothing.

      Every written map now carries the command that made it: ENVI as
      `history = {...}` plus `software = {...}`, TIFF as `ImageDescription`
      and `Software`. Arguments are shell-quoted, and `{`/`}`/newlines are
      replaced rather than escaped -- a mangled history line is better than an
      unparseable header, and `tests/provenance.rs` pins that a path
      containing a brace still leaves the header readable by our own parser.

      No timestamp, deliberately: the same command twice produces identical
      files, which is worth more here than knowing the hour of the run. The
      raster is untouched, so the golden payload comparison is unaffected --
      re-verified byte for byte after the change.
- [x] **12.4 -- `-b`/`-l`/`-B`/`-N`/`-A`.** Half-present is worse than either
      state, so each was taken to one end.

      *Wired up.* `-B band` and `-N low,high` (normality band and interval):
      a region whose centroid in that band falls outside the interval is
      *special* and is held to `Nabsmin` rather than `Nnormin` in Phase 2. The
      logic was already there and correct; there was no way to reach it. Both
      are required together, as in the C, and the C's `high <= 255` becomes a
      check against the input's actual sample range. The two extra auxiliary
      log lines the C prints under `-B` are printed too.

      `-A` allocated the mask, filled it in during Phase 2 and then dropped it
      on the floor. It is now written as `<base>.armask.<pass>`, one uint8
      band, the way `wr_armm` did.

      *Deleted.* `log_band` (`-b`), `lthr` and `lincr` (`-l`) drove the
      per-pass single-band `.log.<n>` debug files -- the same category as
      `-h`, which decision Q6 already dropped. Nothing read the fields. They
      are gone rather than left looking like features.

      `tests/flags.rs` pins the behaviour rather than the plumbing: a bright
      field with 25 isolated dark specks loses them to Phase 2 without
      `-B`/`-N` and keeps them with, and an interval that covers every centroid
      reproduces the no-`-B` run exactly.

---

## 13. Two-input segmentation (Ye et al. 2025)

The 1992 algorithm segments one image. The modification this repo now has to
support segments **two**: Woodcock & Harward's micro-segmentation on Landsat
spectral proxies as before, then a *segment-development* phase that merges those
micro-segments using a **different image over the same grid** — forest structure,
age, or species. Published as Ye, Coops, Wulder & Hermosilla, *ISPRS J.
Photogramm. Remote Sens.* 226 (2025) 381–395; the local copy is in `paper/`,
which is gitignored.

The requirement is that this stays one program. With no second image the
behaviour is exactly what it is today, byte for byte, golden fixtures and all.
Supply a second image and the auxiliary phase is replaced by the
segment-development phase.

### 13.1 What the second phase actually is

Not a re-run of Phase 2 on different pixels — the rules differ:

| | Phase 2 (1992, `seg_apass`) | Segment development (2025) |
|---|---|---|
| Input | the same image as Phase 1 | a second image, different bands, same grid |
| Starting map | the `armap` in progress | the **`rmap`** — normal passes only |
| Merge partner | mutual nearest neighbour | undersized region's nearest neighbour, **made** mutual by write-back |
| Surviving id | the lower id | the **absorbing (smaller)** region's id |
| Distance | f32 | f64 |
| Tie-break | `flip()`, glibc `random()` | none needed (see 13.3) |
| Size rules | `Nabsmin`/`Nnormin`/`Nviable`/`Nmax`/`Nabsmax`, normality band | `Nmin` and `Nmax` only |
| Masking | from the input image's nodata | **from the second image**: a region more than half non-treed is dropped |

The write-back is the substantive idea. §4.2 of the paper: *"for any segment A
that was smaller than the minimum region size, it could be merged with its
nearest neighbour segment B, even if segment A was not the nearest neighbour of
segment B."* Implemented by having A stamp its own distance onto B, so the
mutual-nearest test that follows passes by construction unless a closer A′
claims B first.

### 13.2 The oracle, and the two bugs in it

Elaine's Python (`~/mac2025/segment_python`, commit `427a5a3`) is the definition
of the phase — it produced the published results. It is vendored into
`tools/stage2_oracle/` with a reviewable diff, and `tests/stage2/` holds six
cases generated from it. `tests/STAGE2.md` is the full description; the
generator refuses to write a case it cannot pin.

Two defects had to be dealt with rather than ported:

**The whole-map wipe.** `region_map[tuple(np.array([...]).T.tolist())] = 0` with
an empty list is `region_map[()] = 0` — numpy's spelling of *the entire array*.
Whenever no region was majority-nodata, the phase zeroed every id, found no
adjacencies, merged nothing, and wrote the input back out. It never fired on the
real runs because the structure/age/species layers always contain some
majority-nodata region, so the published results are unaffected; it fires
immediately on clean input, which is why Elaine's own `integ_test.py` was
returning its 5×5 input unchanged and looking like it passed. Guarded in the
vendored copy. **Not** to be reproduced in Rust.

**The unseeded tie-break.** Near-equal candidates are chosen between with
`randint(0, 1)` while iterating a Python `set`, so the phase has no defined
answer. Sweeping all six combinations of (set/ascending/descending order) ×
(keep/take on tie), three of the six fixture cases come out different. Pretending
otherwise would have frozen one arbitrary sample as the oracle. The fixtures pin
a rule instead — **ascending region id, keep the incumbent on a near-tie** —
under which the coin flip is unreachable, so stage 2 needs no RNG at all.

Be exact about the justification, because the 1992 C has two mechanisms that are
easy to conflate. Its *merge survivor* rule is deterministic and id-based:
`if (r < nnbr_id) merge_regions(Spr, r, nnbr_id); else merge_regions(Spr,
nnbr_id, r)`, and `region.c` says so — *"merge the two regions into the region
with the lower REGION_ID"*. Its *candidate tie-break*, in `reg_nnbr`, is
`flip()` — the same unseeded coin the Python uses, commented *"This is biased,
but it does give some randomness to nnbr selection"*. So the chosen rule is
**consistent in spirit with the C's survivor convention, not a reconstruction of
its tie-break**; the C has no deterministic tie-break to reconstruct. The three
tie-insensitive cases would pass under any rule; the three sensitive ones are
what test that this one is implemented.

The oracle enforces it structurally rather than by patching from outside:
`region.py` exposes `ON_TIE` (`keep`/`take`/`random`) and only `random` reaches
`randint`. Regenerating all six cases under the default produced bytes identical
to the committed manifest, which is the proof the flip was never consulted.

A third quirk is faithful and stays: a region's centroid is the mean over **all**
its pixels including nodata ones, so non-treed pixels drag the mean toward zero.
Only the >50 % test excludes them. That is what the published segmentation did.

### 13.3 Specification

Per pass, until a pass merges nothing (region ids are visited in
**raster-first-occurrence order** throughout — the order ids first appear
scanning the region map row by row):

1. Reset every region's recorded nearest distance to +∞.
2. For each region with `npix < Nmin`: scan its bounding box; for each of its
   pixels take the 4-neighbours that are in bounds, unmasked, and a different
   region; among those distinct ids pick the smallest squared Euclidean distance
   between stage-2 centroids, **ascending by id, keeping the incumbent** when
   `|d − dbest| ≤ 1e-6·max(|d|,|dbest|)`. Record `(id, d)`. If `d` is less than
   that neighbour's own recorded distance, write `(this region, d)` onto the
   neighbour.
3. For each region with `npix < Nmin`, in the same order: skip if it or its
   partner already merged this pass, if the partner is 0, if the distance is
   still ∞, if `npix(a) + npix(b) > Nmax`, or if the two recorded distances
   differ by more than `1e-9` relative. Otherwise **a absorbs b, keeping a's id**;
   centroid becomes `(na·ca + nb·cb) / (na + nb)` in f64, in that association.
4. Rewrite the region map.

Before the first pass, centroids are per-band means of the second image over each
region, and every region whose count of all-bands-zero pixels is **strictly more
than half** its area is deleted and its pixels set to region 0.

**On arithmetic.** f64, unlike stage 1's f32 — this is not the C being ported.
The initial means are exactly reproducible without imitating numpy: the samples
are integers, every partial sum is an exact integer below 2^53, so summation
order is irrelevant and one `i64` sum plus one `f64` divide matches bit for bit.

**On adjacency.** The Python maintains a 4-bit-per-pixel adjacency band
incrementally, which forces it to keep a coordinate list per region and makes
merging quadratic. It is equivalent to recomputing the band from the region map
at the start of each pass, because regions participate in at most one merge per
pass so the updated pairs are disjoint. *Verified*: on all six fixtures the
recomputed band is bit-identical to the incremental one at every pass and the
final maps agree. So the Rust keeps no coordinate lists — bbox, `npix`, centroid
and the existing `rband` are enough, and neighbour collection reuses the shape of
`collect_nbrs`/`select_nnbr` already in `segment.rs`.

### 13.4 Shape in this codebase

- `src/stage2.rs`, new. Does **not** extend `Segmenter`: different centroid
  width, different merge direction, no RNG, no contiguity band. Sharing it would
  put a mode flag through the hot loop of a phase that is currently byte-exact.
- `io::envi::read_region_map`, new and separate from `io::envi::read`. Accepts
  ENVI data types 1, 12 and 13 — the fixtures need uint16 (our own writer) and
  uint32 (the NTEMS run). The image reader keeps refusing anything but 8- and
  16-bit; a region map is labels, not DN, and the two should not share a path.
- CLI, additive:

      --stage2 <IMAGE>     second-stage image; enables the phase
      --n2 <NMIN,NMAX>     minimum and maximum region size for it
      --rmap <FILE>        start from this region map, skip stage 1

  Absent `--stage2`, nothing changes. `--rmap` is what lets stage 2 run against
  a region map this program did not produce, which is how four of the six
  fixtures work; it also makes the phase usable on the published NTEMS stage-1
  outputs directly.
- Output stays `<base>.armap.<passes>`, with `passes` counted the Python's way
  (merging passes plus the final no-op one).

### 13.5 Milestones

- [x] **M8a — Test data.** Six cases in `tests/stage2/`, checksum-pinned and
      locked (`tests/verify_stage2.sh`, 43 files), generated from the vendored
      oracle with the arbitrariness sweep as a gate. Coverage across the set:
      every rejection branch (`over_max` 10 715 in `age_capped`, `not_mutual`,
      `busy`, `no_cand`, `inf`), single- and six-band second images, uint16 and
      uint32 region maps, a pre-masked stage-1 map, masking introduced by
      stage 2, and two cases whose stage-1 map this repo's own binary produced.
      *Gate met:* `tests/stage2_fixtures.rs`, 5 tests, deriving the invariants
      from the bytes rather than from the Python — no region split, no invented
      id, masking only grows, nothing over `Nmax` that did not start there.
- [ ] **M8b — `read_region_map` + `--rmap`.** *Gate:* every fixture's region map
      round-trips through the reader at its declared type and count.
- [ ] **M8c — The phase.** *Gate:* all six `expected/armap.<n>` reproduced
      **byte-exactly**, and the pass count matching. Bring up on
      `tiny_synthetic` first — 25 pixels, and it already exercises merge,
      over-`Nmax` and not-mutual rejection.
- [ ] **M8d — Composition.** *Gate:* `--stage2` absent leaves both golden cases
      byte-exact (i.e. `tests/segment_golden.rs` untouched and passing), and a
      single invocation with `--stage2` reproduces `e2e_gsv` from the proxies
      crop without an intermediate file.

### Still open

**`Nviable` has no analogue.** Stage 2 takes only `Nmin`/`Nmax`. The paper's
`Nabsmax` is our `Nmax`; there is no viable-size rule and no normality band. If
you want `-B`/`-N` to apply to the second phase too, say so — it is not in the
Python and I have not invented it.

**Which region should survive a merge.** The Python keeps the *absorbing*, i.e.
smaller, region's id, so ids in the output are not ordered and the map is not
comparable to a 1992 `armap` by id. That is what is implemented, because it is
what the published results were produced with. Changing it to Woodcock's lower-id
convention would invalidate all six expected maps and require regenerating them
from a changed oracle — still possible, but no longer free.
