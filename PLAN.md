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
unambiguous about every decision below. It stays in reserve for §11 Q1.

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

**Consequence for working on a Mac.** Apple's libc `random()` is a *different*
generator. So:

- The Rust implementation must carry its own glibc-compatible `random()` and
  never call the platform's.
- **The C reference built locally will not reproduce the golden files either** —
  not because the port is wrong, but because macOS hands it a different RNG.
  Before trusting the C as an oracle, compile a `glibc_random.c` into it that
  defines `random()`/`srandom()`; the program's own definition wins over libc's
  at link time. Without this, Milestone 0 produces a red herring that looks like
  an algorithm bug.
- Both implementations then share one reference: the same 30 lines, verified
  against a known glibc output vector (seed 1 → `1804289383, 846930886,
  1681692777, 1714636915, 1957747793, …`).

Reproducing the stream also requires the call count and order to match, which
means bit-identical f32 distances (§3.2) and identical neighbour-set order (§3.3).

*Cheap first experiment:* run the shimmed C twice with `flip()` forced to 0 and to
1. If both reproduce the golden bytes, ties never decide anything on these inputs
and the whole problem evaporates. I expect it matters, but it costs ten minutes.

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

- [~] **M0 — C reference on this machine.** (in progress) Build the original; add
      `glibc_random.c` (§3.1); regenerate Case 1 from `misc/temp_byte_bip`.
      *Gate:* the C reproduces `proof/regmap.armap.58` and `regmap.rmap.51`
      byte-exactly. Then run the `flip()`-forced-0/1 experiment and record whether
      the RNG matters at all.
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
- [ ] **M7 — Performance.** Benchmark vs. the C at 250² and 5000²; then parallel
      `reg_nnbr` (§7.3). *Gate:* wall-clock at or below the C on 5000², and M3/M4
      still byte-exact after parallelisation.

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
