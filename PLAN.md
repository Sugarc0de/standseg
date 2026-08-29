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

**Two exits to delete, per your instruction:**
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

`random()` is **never seeded**, so it runs from glibc's default seed of 1. It is
called in `reg_nnbr` whenever a candidate neighbour's distance exactly equals the
running minimum, and (only under `-A`) in `seg_apass`. Early passes over uniform
imagery produce many exact ties, so this is consumed constantly and the output
depends on the whole call sequence.

Reproducing it requires three things to line up:
1. **The generator.** glibc `random()` is TYPE_3: 31-word additive feedback,
   `r[i] = r[i-3] + r[i-31]`, output `(u32)r[i] >> 1`; state seeded by the
   Lehmer recurrence `16807·r[i-1] mod 2³¹-1` (Schrage form), then 310 outputs
   discarded. Roughly 30 lines to port. **But macOS `random()` is a different
   implementation** — if the golden files were produced on a Mac, this port is wrong.
2. **The call count.** Every tie must be detected identically, which means the f32
   distances must be bit-identical (§3.2) and the neighbour set must be in the
   same order (§3.3).
3. **The call order.** Regions are visited in ascending id, serially. Any
   parallelisation must not reorder RNG consumption — see §7.

*Empirical first move, before porting anything:* build the C reference locally and
run it twice with `flip()` hard-wired to 0 and to 1. If both reproduce the golden
bytes, ties never decide anything on these inputs and the RNG is a non-issue.
I expect it does matter, but this is a cheap check that could delete the whole
problem.

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

I would ship the straightforward 9.9 GB version first and treat spilling as a
separate project.

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
  rng.rs         glibc TYPE_3 random() port, behind a trait so BSD can drop in
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

## 9. I/O

**Read:** ENVI (raw + `.hdr`), IPW, TIFF/GeoTIFF, PNG. All must present as
uint8 BIP; the original rejects anything but Byte and the rewrite should too
(§11 Q2). IPW headers are plain text records terminated by a form feed, pixels
byte-aligned at `bytes` per pixel and masked to `bits` — not bit-packed. `nbands`
comes from the `basic_image_i` record, per-band width from `basic_image N`.

**Write:** the region map, sized to the smallest type holding `nreg`
(`nbits ≤ 8` → u8, `≤ 16` → u16, else u32) exactly as `GDAL_write_image` does.
Default output is ENVI raw + `.hdr` so the bytes line up with `proof/`. IPW output
is needed to regenerate the `.armap.58`-style containers.

Geotransform and projection pass through from input to output where the format
carries them (ENVI `map info` / `coordinate system string`, GeoTIFF tags).

---

## 10. Milestones

| # | Deliverable | Gate |
|---|---|---|
| 0 | Build the C reference locally; regenerate Case 1 from `temp_byte_bip`; run the `flip()`-forced-0/1 experiment | Golden bytes reproduced by the C, and we know whether the RNG matters |
| 1 | I/O layer | Read all four formats; round-trip ENVI and IPW byte-exactly |
| 2 | Phase 0 | `nreg` = 55226 (Case 1) and 31609 (Case 2) after `pix_merge` |
| 3 | Phase 1 | Per-pass region counts and merge statistics match `myseg.log` line for line; `rmap.51` and `rmap.17` byte-exact |
| 4 | Phase 2 | `armap.58` and `armap.1` byte-exact. **Definition of done.** |
| 5 | Scale | 15000 × 15000 × 6 completes; memory measured against §6 |
| 6 | Performance | Benchmarked vs. C at 250² and 5000²; parallel `reg_nnbr` (§7.3) still byte-exact |

Milestone 3's gate is the useful one: `myseg.log` records `nreg`, `dmin2`,
`maxpix`, and all seven merge counters for every pass. Diffing those against our
own log localises a divergence to a specific pass and a specific rejection reason,
instead of leaving us with 125000 mismatched bytes and no idea why.

---

## 11. Decisions I need from you

**Q1 — Which machine produced the golden files?**
The `flip()` RNG (§3.1) is glibc-specific. If `.rmap.51` / `.armap.58` came off
the BU Linux cluster, a glibc `random()` port is right. If they were regenerated
on a Mac, I need the BSD implementation instead. The IPW headers record the
command but not the host. *Default if you don't know:* implement glibc, and if
Milestone 3 diverges, try BSD before suspecting the algorithm.

**Q2 — Formats and input types.**
You listed `.tiff`, `.png`, `.ipw`, but the test cases are ENVI and IPW — ENVI is
required whether or not it is on the list. I plan to support ENVI + IPW + TIFF +
PNG natively without linking libgdal. Two sub-questions:
(a) Keep the original's uint8-only restriction, or accept uint16/int16 input by
converting? Case 2's ENVI `_stack` is int16 and the C **rejects** it — that is why
only the `.ipw` works there. Widening input would be new behaviour, not a port.
(b) Should output format mirror the input driver, or always default to ENVI?

**Q3 — Memory ceiling at 15000².**
Peak is ~9.9 GB (§6), and the centroid list cannot shrink without breaking
exactness. Is that acceptable on your target machine? If not, say what ceiling you
need and I will cost out mmap-backed region arrays.

**Q4 — Is byte-exactness required at 15000², or only on the test cases?**
Some attractive optimisations (parallel merge, a proper spatial index instead of
bbox scans, a deterministic tie-break replacing `flip()`) would change output.
If exactness is a test-suite property rather than a program property, I would add
a `--fast` mode and get substantially better large-image performance. Default
assumption: exact always, `--fast` not built.

**Q5 — Keep the 65535-pixel-per-region cap?**
`npix` as `u16` matches the C, keeps the region struct at 12 bytes, and is already
implied by the CLI's `nabsmax ≤ 65535` validation. At 15000² a 65535-pixel region
is 0.03% of the image, so the cap is probably irrelevant — but if you want regions
that can grow arbitrarily, `u32 npix` costs ~1.6 GB more at that size and changes
the `-n` validation rules.

**Q6 — Which flags are in scope for v1?**
The tests exercise only `-t -m -n -o` and 4-way. I plan to also implement
`-8`, `-M` (mask), `-b`/`-l` (log bands), `-B`/`-N` (normality), and `-A`
(auxiliary mask) since they are woven through the algorithm and cheap to carry.
I plan to **drop** `-h` (writes `.cband`/`.rlist` files for `hsegment`, a program
that is not in this repo) and `-S` (per your instruction). Confirm, particularly
on `-h`.
