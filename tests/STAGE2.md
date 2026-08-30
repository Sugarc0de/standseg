# Stage-2 fixtures — the two-input segment-development phase

`tests/stage2/` is the oracle for the **modified** algorithm: Ye et al. 2025
(ISPRS J. Photogramm. Remote Sens. 226, 381–395), which keeps Woodcock &
Harward's micro-segmentation as stage 1 but replaces the auxiliary phase with a
*segment-development* phase driven by a **second, different image** — forest
structure, age, or species — over the same grid.

These fixtures are **not** `tests/golden/`. Read the distinction before touching
either:

| | `tests/golden/` | `tests/stage2/` |
|---|---|---|
| What it is | 1992 reference outputs, captured from the original program | Outputs of Elaine's Python stage-2 implementation, regenerated here |
| Status | **Inviolable.** Never regenerate. | Regenerable from `tools/stage2_oracle/`, deliberately and with the diff reviewed |
| Pins | The port is faithful to the C | The Rust stage 2 is faithful to the Python |

Verify at any time: `tests/verify_stage2.sh`

## Where the data comes from

Tile 399, Alberta, from the NTEMS 2019 update — the same product family the paper
used (`/Volumes/easystore/UBC/first_project/ntems_2019_update/ab/processed_tiles/tile_399`,
an external drive, so the fixtures are the durable copy, not the tile).

- **Stage-1 region maps.** Either a crop of the real
  `proxies_t_50_m_0.2_n_9_18_36_tile_399.rmap.26` (the published run: `-t 50
  -m 0.2 -n 9,18,36` on the 6-band Landsat BAP proxies), or — for the `e2e_*`
  cases — produced here by **this repo's own Rust binary** on a crop of the same
  proxies. The `e2e_*` cases are what pin that the two stages compose.
- **Stage-2 layers.** `elev_p95`, `elev_cv`, `gross_stem_volume` (structure),
  `age`, and the 6-band `species` probability stack. All are the `-norm` variants:
  uint8, rescaled to 1–255 with **0 reserved for non-treed**, exactly as §4.2 of
  the paper describes.

Note that stage 1 consumes the **`rmap`** (normal passes), not the `armap` — the
whole point is that stage 2 replaces the auxiliary phase.

## The six cases

| case | grid | stage-2 layer | Nmin | Nmax | what it is for |
|---|---|---|---|---|---|
| `tiny_synthetic` | 5×5×1 | synthetic | 4 | 9 | Elaine's own `integ_test.py` example. Hand-checkable; exercises merge, over-Nmax and not-mutual rejection in 25 pixels. |
| `p95_250` | 250×250×1 | `elev_p95` | 80 | 8000 | The paper's parameters on a single structural band. 31 % non-treed. |
| `species_250` | 250×250×6 | `species` | 80 | 8000 | Multi-band. Zero is a legitimate probability here, so it pins the *all-bands-zero* reading of the invalid-pixel rule. |
| `age_capped` | 250×250×1 | `age` | 60 | 200 | Deliberately tight `Nmax` — the only case where the maximum-region-size rejection actually binds (10 715 times). |
| `e2e_gsv` | 200×200×1 | `gross_stem_volume` | 50 | 8000 | Stage-1 map produced by our own Rust. uint16 region map. |
| `e2e_masked` | 200×200×1 | `elev_cv` | 50 | 8000 | Stage 1 run under `-M`, so the input region map **already contains masked pixels** before stage 2 adds its own. |

Each case directory holds `input/rmap`, `input/layer`, `expected/armap.<passes>`
(all ENVI, BSQ, little-endian, `.hdr` alongside) and a `case.json` with the
parameters, the crop window, region counts, and per-pass merge/rejection
counters. `cases.json` at the top is all six.

Region maps are uint32 (ENVI data type 13) for the real-tile crops and uint16
(type 12) for the `e2e_*` cases, because that is what each producer emitted —
which means **the region-map reader has to accept 1, 12 and 13**. That is a
different code path from the image reader, which deliberately still refuses
anything but 8- and 16-bit.

## The algorithm these fixtures pin

Per pass, until a pass merges nothing:

1. Every region's recorded nearest-neighbour distance is reset to +∞.
2. For each region with `size < Nmin`, in **raster-first-occurrence order** of
   region ids: scan its bounding box for pixels with an adjacent pixel of a
   different, unmasked region; among those neighbours pick the smallest squared
   Euclidean distance between **stage-2 centroids**. Record it, and if that
   distance beats the neighbour's own recorded distance, **write it back onto the
   neighbour**. That write-back is the paper's relaxation: an undersized segment
   may merge with its nearest neighbour even when it is not that neighbour's
   nearest.
3. In the same order, merge: skip if either side is already involved this pass,
   if the partner is 0, if the distance is still ∞, if `size(a) + size(b) > Nmax`,
   or if the two recorded distances no longer match. Otherwise **the small region
   absorbs the larger one and keeps its own id** — the opposite of stage 1, which
   merges into the lower id.
4. Centroids combine size-weighted; the region map is rewritten at end of pass.

Before any of that, centroids are means of the **stage-2** samples over each
region, and any region whose all-bands-zero pixel count is **strictly more than
half** its area is dropped and its pixels set to region 0. That is how non-treed
area enters: `> 50 % non-treed ⇒ excluded`, §4.2 of the paper.

### Arithmetic

Centroids and distances are **f64** here, unlike stage 1's f32 — this is a
different algorithm, not the C being ported, and the Python that defines it uses
Python floats throughout.

The initial means are reproducible without heroics: numpy sums uint8 samples in a
float64 accumulator, every partial sum is an exact integer below 2^53, so
pairwise and sequential summation give the same value and a single correctly
rounded division finishes it. An `i64` sum and one `f64` divide in Rust matches
bit for bit. `update_centroids` is then
`(n1*c1 + n2*c2) / (n1+n2)` in f64, in that association.

Tie comparison is `math.isclose(a, b, rel_tol=1e-6)`, i.e.
`|a-b| <= 1e-6 * max(|a|,|b|)`; the mutual-distance check uses the default
`rel_tol=1e-9`.

### The tie-break — read this before implementing

The Python picks among near-equal candidates with an unseeded `randint(0, 1)`,
iterating a Python `set` whose order is an implementation detail. Two arbitrary
choices, neither reproducible in Rust.

The generator therefore sweeps all six combinations (`set`/`asc`/`desc` order ×
keep/take on tie) for every case and records whether they agree. They do **not**
always agree: `age_capped`, `e2e_gsv` and `e2e_masked` are sensitive to it. So
the fixtures pin a rule rather than pretending the question does not exist:

> **Visit candidate neighbours in ascending region id; on a near-tie keep the
> incumbent.** The smallest id among near-equal candidates wins.

Under this rule the coin flip is never reached — **stage 2 needs no RNG at
all**, unlike stage 1, which needs glibc `random()` call for call. Proved rather
than assumed: `tools/stage2_oracle/region.py` now selects the tie policy with
`ON_TIE` (`keep`/`take`/`random`) and only the `random` setting touches
`randint`, and regenerating all six cases with the default produced files
byte-identical to the ones already in `STAGE2.sha256`.

**This is a decision, not a reconstruction, and it is worth being exact about
what the 1992 C does — because two different mechanisms are easy to conflate:**

| the C | mechanism |
|---|---|
| which region *survives a merge* | the **lower id**, deterministic: `if (r < nnbr_id) merge_regions(Spr, r, nnbr_id); else merge_regions(Spr, nnbr_id, r);` — `region.c` even says so in the comment, *"merge the two regions into the region with the lower REGION_ID"* |
| which of several *equidistant candidates* to pick | `flip()`, i.e. unseeded `random() & 01` — `reg_nnbr` in `region.c`, commented *"This is biased, but it does give some randomness to nnbr selection"* |

So the C's tie-break is **not** a lower-id rule; it is a coin flip, exactly like
the Python's. Ascending-id-keep-incumbent is chosen because it is deterministic,
cheap in both languages, and consistent in spirit with the C's *survivor* rule —
not because either implementation did it. The three tie-insensitive cases would
pass under any rule; the three sensitive ones are what actually test that this
one is implemented.

## Regenerating

`tools/stage2_oracle/` holds the exact Python used, vendored from
`~/mac2025/segment_python` at commit `427a5a3`. `diff` it against upstream: the
changes are instrumentation, an `O(1)` membership set in place of an `O(n)` list
scan in `update_adjacent_regions` (pure optimisation — the function only ORs
bits, so neither coords order nor lookup structure can change its result; the
generator asserts both spellings agree), the order/flip switches above, and one
guarded bug:

    if majority_invalid_region_ids:        # ADDED
        region_map[tuple(np.array([...]).T.tolist())] = 0

With an empty list that index expression is `region_map[()]` — **the whole
array** — so the stage silently wiped every region id and did nothing whenever no
region was majority-nodata. It never fired on the real runs, because the
structure/age/species layers always contain some majority-nodata region; it fires
immediately on clean synthetic input, which is why `tiny_synthetic` was returning
its input unchanged. See PLAN.md §13.2.

Regeneration needs the external drive mounted and a venv with
`numpy rasterio numba`:

    python tools/stage2_oracle/gen_fixtures.py     # writes tests/stage2/
    python tools/stage2_oracle/verify_fixtures.py  # invariant check, independent code

## Comparing a rewrite run

Same discipline as `tests/golden/`: **program output goes to `build/out/`**,
never into `tests/stage2/`.

    fast_segment --rmap tests/stage2/<case>/input/rmap \
        --stage2 tests/stage2/<case>/input/layer --n2 <Nmin>,<Nmax> \
        -o <case> --outdir build/out
    cmp build/out/<case>.armap.<n> tests/stage2/<case>/expected/armap.<n>

`cargo test` does this for all six, through the library
(`tests/stage2_match.rs`) and through the command line
(`tests/stage2_cli.rs`).

Byte-exact equality of the region-map payload is the pass condition. Region ids
carry meaning here — the surviving id is the *absorbing* region's — so "same
partition, different numbering" is a failure, not a near miss.

`case.json` also records the oracle's per-pass `considered` / `no_cand` / `busy`
/ `inf` / `over_max` / `not_mutual` / `merged` counts, and those are checked too.
They are what localises a divergence to a pass and a reason instead of a byte
offset — and two ways of getting the map right while getting them wrong are
written up in PLAN.md §13.7.
