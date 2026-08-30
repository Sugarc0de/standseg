# fast_segment

A Rust reimplementation of Harward & Woodcock's nested-hierarchical
region-growing image segmenter — the `segment` program used for tree stand
delineation from satellite imagery.

It reproduces the original C program's output **byte for byte**, handles images
the C cannot (it segfaults above roughly 5000 × 5000), and runs about 2× faster
than a `-O2` build of the original.

> Woodcock, C. and V. J. Harward. 1992. *Nested-hierarchical scene models and
> image segmentation.* International Journal of Remote Sensing 13(16): 3167–3187.

## Quick start

Needs only a Rust toolchain — no GDAL, no system libraries.

```bash
git clone https://github.com/Sugarc0de/fast_segment
cd fast_segment
cargo build --release
```

Segment the bundled 250 × 250 four-band test scene:

```bash
mkdir -p out
./target/release/fast_segment \
    -t 10 -m .1 -n 15,15,100,2500,2500 \
    -o demo --outdir out \
    tests/golden/misc/temp_byte_bip
```

That writes `out/demo.rmap.51` and `out/demo.armap.58` (plus `.hdr` sidecars).
Those bytes are identical to the reference output the original C produced,
which ships in this repo:

```bash
cmp out/demo.rmap.51  tests/golden/test_3456/expected/proof/regmap.rmap.51
cmp out/demo.armap.58 tests/golden/test_3456/expected/proof/regmap.armap.58
```

Or just run the test suite, which checks that and a good deal more:

```bash
cargo test --release
```

## What it does

Two phases, both in a single run.

**Phase 1 — normal passes.** Every region finds its nearest neighbour in
spectral space (Euclidean distance between band centroids). Mutually-nearest
pairs within the tolerance merge, one merge per region per pass, until a pass
produces none. Writes `<base>.rmap.<pass>`.

**Phase 2 — auxiliary passes.** Regions still below the minimum size are forced
to merge with their nearest neighbour, this time with no distance ceiling, until
no region can. Writes `<base>.armap.<pass>`.

The pass number is part of the filename, so `demo.armap.58` means the auxiliary
phase converged after 58 passes.

### A second phase driven by a second image

Woodcock & Harward segment one image. Ye et al. (2025) keep phase 1 as a
*micro-segmentation* over Landsat spectral proxies, then replace phase 2 with a
*segment-development* phase that merges those micro-segments using a **different
image over the same grid** — forest structure, age or species. That variant is an
option here, not a fork: with no second image the program does exactly what it
did before, byte for byte.

It stays one command. Add `--stage2` and the same run micro-segments the
proxies, then develops the result against the second image:

```bash
# one image -- Woodcock & Harward, unchanged
fast_segment -t 50 -m 0.2 -n 9,18,36 -o stands proxies

# two images -- the same run, plus the second phase
fast_segment -t 50 -m 0.2 -n 9,18,36 \
    --stage2 elev_p95 --n2 80,8000 -o stands proxies
```

The second writes both maps: `stands.rmap.41` is the stage-1 micro-segmentation,
kept exactly as the one-image run would have written it, and `stands.armap.39` is
the developed result. The intermediate is a real output, not a temporary — you
can hand it back later:

```bash
# already have a stage-1 map? develop it against a different second image
fast_segment --rmap stands.rmap.41 --stage2 age --n2 60,8000 -o stands_age
```

| Option | Meaning |
|---|---|
| `--stage2 <image>` | The second image. Same grid, different data. Enables the phase, replacing the auxiliary one. |
| `--n2 <Nmin,Nmax>` | Size rules for it: merge regions up to `Nmin`, never across `Nmax`. |
| `--rmap <file>` | Optional shortcut: take stage 1's region map from a file and skip stage 1. Then no input image or `-t` is needed. |

Two things differ from phase 2 that are easy to trip over. A region merges with
its nearest neighbour **even when it is not that neighbour's nearest** — the
paper's relaxation, and the reason small segments get absorbed at all. And the
surviving id is the *small* region's, not the lower one's, so ids in a
`--stage2` map are not comparable to a 1992 `armap` by id.

A region whose second-stage pixels are **more than half zero** is dropped and its
pixels set to 0. That is how non-treed area enters: the layers are rescaled to
1–255 with 0 reserved for non-treed, so the mask comes from the second image, not
from `-M`.

The oracle is Elaine Ye's Python implementation, the one behind the published
results. Six cases generated from it live in `tests/stage2/`, and all six
reproduce byte for byte, at the same pass counts and with the same per-pass merge
and rejection counts. `tests/STAGE2.md` describes them and the two bugs in that
Python that had to be dealt with rather than ported; `PLAN.md` section 13 has the
design.

## Usage

```
fast_segment [OPTIONS] -t <TOLS> -o <BASE> <IMAGE>
```

| Option | Meaning |
|---|---|
| `-t <t1,t2,…>` | Segmentation tolerances. Required. Each produces its own region map. |
| `-o <base>` | Basename for output files. Required. |
| `-m <cm>` | Merge coefficient, `0 < cm ≤ 1`. Caps merges per pass at `cm × nregions`. Default 1 (no cap). |
| `-n <a,b,c,d,e>` | `Nabsmin,Nnormin,Nviable,Nmax,Nabsmax` — region size rules, see below. |
| `-8` | Use 8-way connectivity instead of 4-way. |
| `-M <file>` | Mask image; pixels valued 0 are excluded. |
| `-B <band>` | Zero-based band carrying the normality criterion. Requires `-N`. |
| `-N <low,high>` | Normality interval. A region whose `-B` band centroid falls outside it is *special*, and Phase 2 holds it to `Nabsmin` instead of `Nnormin` — the way to stop small non-forest patches being absorbed into the stands around them. |
| `-A` | Also write the auxiliary region map mask, `<base>.armask.<pass>`: 0 where Phase 2 absorbed a region, 1 elsewhere. |
| `--nodata <v>` | Treat pixels with this value as nodata. May be negative (Landsat's `-9999` fill). |
| `--nodata-any` | A pixel is nodata if *any* band matches, rather than all bands. |
| `--outdir <dir>` | Where to write output. Default `.`. |
| `--format <envi\|tiff>` | Output format. Default `envi`. |
| `--threads <n>` | Worker threads. `0` = one per core, `1` = serial. Output is identical either way. |

### The `-n` size rules

```
-n Nabsmin,Nnormin,Nviable,Nmax,Nabsmax
```

- **Nabsmin** — floor for "special" regions (those outside the `-B`/`-N` normality interval).
- **Nnormin** — floor for ordinary regions. Phase 2 exists to enforce this; if it is 1, Phase 2 has nothing to do.
- **Nviable** — two regions may not merge if *both* already have this many pixels.
- **Nmax** — size ceiling during Phase 1.
- **Nabsmax** — absolute ceiling, relaxed to during Phase 2.

They must satisfy `0 < Nabsmin ≤ Nnormin ≤ Nviable ≤ Nmax ≤ Nabsmax`, and `0`
in any position means "no limit".

The original stopped at 65535 in every position, because its pixel counter was
an `unsigned short` — so a default run silently stopped growing a stand at
65535 pixels, which at 1 m resolution is a 256 m square. That ceiling is gone:
values above 65535 are accepted and "no limit" now means no limit. Ask for
`-n ...,65535,65535` if you want the old behaviour back.

## Image formats

Read: **ENVI** (raw + `.hdr`), **IPW**, **TIFF/GeoTIFF**, **PNG**.
Write: **ENVI** (default, byte-compatible with the original) or **TIFF**.

Format is detected from content first, extension second — the reference input
`temp_byte_bip` has no extension at all.

**Sample widths: 8-bit unsigned, 16-bit unsigned, 16-bit signed.** The 1992
original was uint8-only, which is what an 8-bit TM scene was. Landsat 8/9 and
Sentinel-2 are 12-bit data delivered in a 16-bit container (int16 for
Collection 2 surface reflectance), and rescaling that to a byte before
segmenting discards radiometry that changes where the boundaries land — so the
wide types are read as they are. 32-bit and floating-point samples are still
refused: tolerances and distances here are integer DN.

Widening is a generalisation, not a second algorithm. Feeding the same values
in as `u8`, `u16` or `i16` produces bit-identical region maps
(`tests/wide_input.rs`), and the 8-bit path is unchanged — the golden fixtures
still reproduce byte for byte.

One consequence worth stating: **tolerance is in DN**, so it does not carry
across widths. `-t 10` on 8-bit reflectance is roughly `-t 350` on the same
scene at 16 bits. There is no automatic scaling; pick the tolerance for the
data you have.

Bands map to samples-per-pixel, so an RGB image is 3 bands and a 6-band
satellite stack is 6 bands. For PNG, note that an alpha channel reads as an
ordinary band and will take part in the spectral distance — use `--nodata` or
`-M` for transparency instead.

### Provenance

Every map records the command that produced it — ENVI as `history` and
`software` keys, TIFF as the `ImageDescription` and `Software` tags. This is
what IPW did in 1992, and it is the only reason the invocation behind the
reference outputs was still recoverable eleven years later:

```
history = {fast_segment -t 10 -m .1 -n 15,15,100,2500,2500 -o stands scene.tif}
software = {fast_segment 0.1.0}
```

There is no timestamp, on purpose: running the same command twice produces
byte-identical files.

## Nodata (water, non-treed area)

Nodata pixels are excluded from segmentation entirely: they get region 0, never
join a region, never contribute to a centroid, and no region may grow across
them. This matters for stand delineation — without it, stands merge across
lakes.

Nodata comes from, in order of precedence:

1. `--nodata <value>` on the command line;
2. the file's own declaration (ENVI `data ignore value`, GeoTIFF `GDAL_NODATA`);
3. an explicit `-M` mask, which is combined with either of the above.

By default a pixel is nodata when **all** bands equal the value — masked-to-land
imagery carries 0 across every band over water, while a legitimate 0 in a single
band is ordinary dark ground that should still segment. `--nodata-any` switches
to the other reading.

Output headers declare `data ignore value = 0` so nodata round-trips.

```bash
fast_segment -t 10 -m .1 -n 15,15,100,2500,2500 \
    --nodata 0 -o stands --outdir out scene.tif
```

## Performance

Measured on an M-series laptop, 10 cores, 6-band imagery:

| scene | C `-O0` (original Makefile) | C `-O2` | this, serial | this, parallel |
|---|---|---|---|---|
| 5000 × 5000 | 24.5 s | 12.6 s | 13.3 s | **9.7 s** |
| 15000 × 15000 | *segfaults* | *segfaults* | 157.9 s | **77.6 s** |

A 5000 × 5000 tile needs about 0.7 GB; 15000 × 15000 peaks around 5 GB. Every
run prints its own array-size breakdown before the first pass.

The C fails above roughly 5000 × 5000 because `ecalloc` takes a signed 32-bit
byte count and the centroid list overflows it. This version uses `usize`
throughout.

Only the nearest-neighbour scan is parallel; the merge loop is inherently
sequential, so expect about 2× rather than one-per-core. That is the price of
bit-exact output, and it is deliberate.

## Is it really identical?

Fidelity rests on three independent checks, all in `cargo test`:

1. **The golden fixtures.** Both reference cases reproduce byte for byte, both
   phases, converging on the same pass numbers (51/58 and 17/1). Every
   per-pass statistic in the original's log matches.
2. **Byte-identity with the C at 5000 × 5000** — 400× the area of the test
   cases, on data neither implementation was tuned for.
3. **Serial and parallel agreement**, with the golden cases forced through the
   parallel path.

Getting there required reproducing some details that are easy to miss: the
program breaks distance ties with an unseeded `random()`, and that tie-break is
load-bearing — forcing it to a constant changes the output. See `PLAN.md §3` for
the full list of hazards.

## Repository layout

```
src/                  the segmenter
tests/golden/         1992 reference inputs and outputs, checksum-pinned (read-only)
tests/stage2/         two-image segment-development fixtures, checksum-pinned
tests/                integration tests, incl. both byte comparisons
tools/stage2_oracle/  the Python that defines the second phase, vendored
reference/csegment/   the original C, buildable as a debugging oracle
PLAN.md               design notes: algorithm, port hazards, memory, milestones
```

`tests/golden/` is the oracle and is pinned by `tests/GOLDEN.sha256`. Git does
not preserve read-only permissions, so after cloning you may want:

```bash
tests/lock_golden.sh      # make the fixtures read-only
tests/verify_golden.sh    # confirm nothing has drifted
```

`tests/stage2/` works the same way (`tests/lock_stage2.sh`,
`tests/verify_stage2.sh`) with one difference: it is *regenerable*, from
`tools/stage2_oracle/`, whereas `tests/golden/` never is.

The C in `reference/csegment/` is not needed to build or use this program. It
exists so a future divergence can be debugged by instrumenting the original.
See `reference/csegment/PORTING.md` for what had to change to build it on macOS
— notably a genuine undefined-behaviour bug in the original's `set.c`.

## Credits

Original algorithm and C implementation by Jud Harward and Curtis Woodcock,
Boston University, building on the IPW library from UC Santa Barbara.
