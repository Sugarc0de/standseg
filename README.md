# fast_segment

[![CI](https://github.com/Sugarc0de/fast_segment/actions/workflows/ci.yml/badge.svg)](https://github.com/Sugarc0de/fast_segment/actions/workflows/ci.yml)

Multi-resolution region-growing segmentation for raster imagery, in Rust.

It implements the algorithm from Ye et al. (2025), which extends Harward and
Woodcock's 1992 segmenter to take a **second data layer**. The original grows
regions in one image and then forces the undersized ones to merge with whatever
is nearest. The modification keeps that first stage as a micro-segmentation, and
replaces the second with one that develops those micro-segments against a
different image over the same grid — forest structure, age, species, anything on
the same pixels.

The original algorithm is still in here and still exact. With one image the
program does what the 1992 C did, byte for byte. The second layer is an option,
not a fork, so both papers' methods run from the same binary.

It also handles images the C cannot (it segfaults above roughly 5000 × 5000) and
runs about 2× faster than a `-O2` build of it.

It was built for forest stand delineation from Landsat, and that is what it has
been tested on. Nothing in the algorithm is forest-specific, though. The
exclusion rule is just "drop any region that is more than half nodata"; in our
data those zeros were non-treed area, but they can be cloud, water, or any mask
you like.

> Ye, E., N. C. Coops, M. A. Wulder and T. Hermosilla. 2025. *A multi-resolution
> forest stand segmentation algorithm integrating Landsat imagery and forest
> structural, age, and species attributes.* ISPRS Journal of Photogrammetry and
> Remote Sensing. https://doi.org/10.1016/j.isprsjprs.2025.05.023
>
> Woodcock, C. and V. J. Harward. 1992. *Nested-hierarchical scene models and
> image segmentation.* International Journal of Remote Sensing 13(16): 3167–3187.

## Why I wrote this

I did my master's at the Faculty of Forestry at UBC, in the remote sensing lab.
The work needed forest stands delineated from Landsat, and the tool for that in
our lab was `segment`. The algorithm is from 1992 and it is still good. Three
things stood in the way of using it.

The first is size. The C segfaults somewhere above 5000 × 5000 pixels, because
`ecalloc` takes a signed 32-bit byte count and the centroid list overflows it.
Our tiles are 5000 × 5000, so I was already at the edge, and the national
coverage is thousands of them.

The second is that it no longer builds cleanly. Getting it to compile on macOS
turned up a genuine undefined-behaviour bug in its own `set.c`. You can work
through that yourself, but it is not something to hand a collaborator who just
wants to segment an image.

The third is about my own work rather than the C. The point of my paper was that
one image is not enough to find a stand. Spectral response tells you where the
canopy changes; it does not tell you that two patches with the same reflectance
are different ages or different species. So I changed the second stage to merge
micro-segments using a second layer — structure, age, species — instead of just
absorbing whatever was smallest. That is the contribution, and the C has no way
to do it.

I wrote that version in Python. It gives the right answer and takes about 25
minutes and 6 GB per tile. That is fine for a thesis. It is not fine for
something meant to run over a country, and it meant the method existed but was
not really usable by anyone else.

There is a fourth reason I only understood while doing the rewrite. Both programs
decide near-ties with a coin flip, and ties are not rare. On a single-band 8-bit
layer, which is what most of my runs used, close to half of all pixels have more
than one nearest neighbour in the first phase (`PLAN.md` §13.7 has the numbers).

The two coins are not equally bad. The C's is `random() & 01`, never seeded, so it
runs from the default seed and gives the same sequence every time. The answer is
arbitrary, but it repeats, and it repeats across machines here because this
version ports glibc's generator instead of calling the platform's.

My Python is the real problem. It called `randint(0, 1)` while iterating a Python
`set`, so the second phase has no defined answer at all. Sweeping the six
reasonable ways to resolve it, three of my six test cases come out different. So
the second-phase maps from my 2023 runs cannot be regenerated from their own
inputs, by me or by anyone else, and I did not know that when I published. This
version pins the rule — ascending region id, keep the incumbent on a near-tie —
and under that rule the coin is never reached, so the second phase needs no
random numbers at all.

So the reasons are, in order: it has to handle a real tile, it has to build, it
has to be fast enough for a national run, and it has to give the same answer
twice.

## Who this is for

If you read the 2025 paper and want to run the method, this is it. The Python
behind the paper is vendored here as the reference the Rust is checked against,
but this is the version to actually use.

More generally, people doing stand delineation or similar segmentation on
satellite imagery who want something free, scriptable, and checkable. Of course there is commercial
software for this, and eCognition is what most people in remote sensing reach
for. It is a good tool. It is also expensive, closed, and there is no way to
check its output against a reference, which is the part that matters if you are
publishing the result.

I am not claiming this segments better than eCognition. It uses a different
algorithm and I have not compared them. What I am claiming is that you can read
it, run it for free, batch it, and verify it — and that the second phase is
something eCognition does not do at all.

## What this is not

There is no GUI. It does not classify anything, compute object features, or
export a multi-scale hierarchy. It reads a raster and writes a region map. If you
want to draw polygons by hand or tune a scale parameter with a slider, this is
the wrong tool and QGIS or eCognition is the right one.

## Install

**A prebuilt binary.** Grab the one for your platform from the
[releases page](https://github.com/Sugarc0de/fast_segment/releases), unpack it,
and run it. There is nothing to install: one static-ish executable, no GDAL, no
Python, no system libraries.

**With Rust, from source.**

```bash
cargo install --git https://github.com/Sugarc0de/fast_segment
```

**Or clone it,** which is what you want if you also intend to run the tests:

```bash
git clone https://github.com/Sugarc0de/fast_segment
cd fast_segment
cargo build --release
```

## Quick start

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

The second image may be **32-bit float**, which is how structural layers
(canopy height, biomass, age, z-scores) normally ship. The first image may not:
stage 1's distances and tolerances are integer DN, and that restriction is
enforced by the type system rather than a check that could be forgotten. Pass a
float raster as the input and the program says so and points at `--stage2`.

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

## Choosing `-t`

`-t` is a distance in **raw DN units**, not a percentage and not a scale factor.
This is the one parameter that will waste your afternoon, because getting it
wrong does not produce an error. It produces a map.

The published parameters (`-t 50 -m 0.2`) are for layers rescaled to 0–255.
Modern imagery is not: Landsat Collection 2 and Sentinel-2 both ship as 16-bit,
where the same number means something roughly 250× smaller. On a 16-bit Landsat
stack running 0–8990:

| tolerance | what the first pass does | regions in the final map |
|---|---|---|
| `-t 10` | merges almost nothing — 62 498 of 62 500 pixels stay separate | 5 304 |
| `-t 350` | merges normally — 55 056 of 62 500 | 3 983 |

At `-t 10` the spectral phase is inert, the size rules force merges anyway, and
you get a plausible-looking map that is 33 % off and shaped by region size rather
than by the image. The program now notices this and warns on stderr, but the
warning is a safety net, not a substitute for scaling the number.

A reasonable starting point is to scale your tolerance by the ratio of the data
ranges. Going from 8-bit to a 0–8990 stack is a factor of about 35, which is how
`-t 10` becomes `-t 350` above.

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
wide types are read as they are. Wider integers and 64-bit float are still
refused.

**32-bit float, for the second stage only.** Structural layers ship as float, so
`--stage2` reads them; stage 1 does not, because its tolerances and distances are
integer DN. Reproducing the reference implementation on a float layer means
reproducing *numpy's* float32 accumulation, including the fact that it sums a
one-band region contiguously (pairwise) and a multi-band one strided
(sequentially). Summing in f64 instead — which is more accurate — moves 5.7 % of
the output pixels. Both orders are implemented and checked against numpy over 285
cases; see PLAN.md §13.8.

Widening is a generalisation, not a second algorithm. Feeding the same values
in as `u8`, `u16` or `i16` produces bit-identical region maps
(`tests/wide_input.rs`), and the 8-bit path is unchanged — the golden fixtures
still reproduce byte for byte.

Because tolerance is in DN, it does not carry across widths — see
[Choosing `-t`](#choosing--t) above.

Bands map to samples-per-pixel, so an RGB image is 3 bands and a 6-band
satellite stack is 6 bands. For PNG, note that an alpha channel reads as an
ordinary band and will take part in the spectral distance — use `--nodata` or
`-M` for transparency instead.

### What is not read, and why

**Satellite product packages.** A Sentinel-2 `.SAFE` directory, a Landsat tar
bundle, a `.jp2` — none of these are read, and that is a scope decision rather
than a gap. This program takes *a raster of values on a grid*: spectral bands, or
a structural attribute layer. Unpacking vendor product formats, applying scale
factors and offsets, resampling 10/20/60 m bands onto a common grid, and reading
cloud masks are all jobs GDAL already does well, and doing them badly here would
be worse than not doing them. Convert first, then segment:

```bash
gdal_translate B04.jp2 B04.tif        # then feed the .tif to fast_segment
```

**64-bit float and 32-bit integer samples.** Not implemented. Both are refused
with a message naming the type rather than failing obscurely. Nothing in
practice needs them: Landsat and Sentinel-2 are uint16, HLS and MODIS are int16,
structural layers are float32, and a header survey of 187 archived experiment
inputs found 186 uint8 and one float32. If you hit one, it is a small addition —
open an issue.

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

## Nodata (water, cloud, non-treed area)

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

### The two-image variant, on a whole NTEMS tile

Tile 399 (Alberta) at its native size: 5000 × 5000 × 6 Landsat proxies for
stage 1, 5000 × 5000 × 1 `elev_p95` for stage 2, with the paper's parameters
(`-t 50 -m 0.2 -n 9,18,36`, `--n2 80,8000`).

| run | wall | peak RSS |
|---|---|---|
| one image, both original phases | 31.0 s | 1.17 GB |
| **one command, `--stage2`** | **35.6 s** | **1.20 GB** |
| segment development alone, from a saved `.rmap` | 11.5 s | 0.53 GB |

Stage 1 leaves **4 377 977** micro-segments; segment development takes them to
**122 552** stands in 114 passes, excluding 1 863 704 regions that were more than
half non-treed. It adds about 4 s and 30 MB to a run that was already happening —
the phase never holds the image, only centroids and bounding boxes — and the
`.rmap.69` it starts from is byte-identical whether or not `--stage2` was asked
for.

Against the Python that defines the phase, on a 1000 × 1000 crop (167 293
input regions): **0.28 s versus 14.5 s**, same 77 passes, byte-identical output.

## How this was checked

A segmentation is hard to eyeball. Two maps can look the same and be different
everywhere it matters, so "close enough" was not allowed to count here. The rule
for the whole rewrite was byte equality against a program that already worked,
and there were two of those: the original C for the first phase, and my Python
for the second.

**The original algorithm, against the C.** Both 1992 reference cases reproduce
byte for byte, in both the normal and auxiliary phases, at the same pass numbers
(51/58 and 17/1), and every per-pass statistic in the original's log matches. The fixtures ship in `tests/golden/`,
checksum-pinned and read-only, so you can check this yourself:

```bash
cargo test --release
```

There is also a byte-identical run against the C at 5000 × 5000, which is 400×
the area of the test cases, on data neither implementation was tuned for.

**Segment development, against the Python.** That phase has no 1992 oracle, so the
Python that produced the published results is the reference. Six generated cases
in `tests/stage2/` reproduce byte for byte, at the same pass counts and the same
per-pass merge and rejection counts. Beyond those, four whole 5000 × 5000 NTEMS
tiles were run through both implementations and compared with `cmp`:

| tile | second-stage layer | passes | result |
|---|---|---|---|
| 397 | 6-band species, uint8 | 134 | identical, 100 000 000 bytes |
| 473 | 3-band structure, uint8 | 154 | identical |
| 474 | 1-band biomass, uint8 | 139 | identical |
| 219 | 3-band age, float32 | 130 | identical |

That is 400 MB of region map with no differing bytes. The Rust runs took 21–26 s
each; the Python took 1440–1780 s and 4.3–7.2 GB.

**What does not reproduce, and why you should know.** Neither implementation
reproduces my 2023 second-phase maps, and neither can, because of the undefined
tie-break described above. The two disagree with those maps on the *same* pixels
— about 1 % on single-band layers, 0.02 % on six-band ones — which is what a
tie-driven difference looks like and is how I know it is the tie-break rather
than a bug.

One part does reproduce exactly: the set of dropped nodata regions has no
tie-break in it, and it matches the 2023 maps on all 52 experiments tested, 0
differing pixels.

If you are using this for published work, that is the honest summary: the method
reproduces, the specific 2023 maps do not, and from here on a run is repeatable
because the tie rule is fixed and written down.

Serial and parallel paths are also checked against each other, with the golden
cases forced through the parallel path. `PLAN.md` §3 lists the details that were
easy to miss.

## Repository layout

```
src/                  the segmenter
tests/golden/         1992 reference inputs and outputs, checksum-pinned (read-only)
tests/stage2/         two-image segment-development fixtures, checksum-pinned
tests/                integration tests, incl. both byte comparisons
tools/stage2_oracle/  the Python that defines the second phase, vendored
reference/csegment/   the original C, buildable as a debugging oracle
PLAN.md               design notes: algorithm, port hazards, memory, milestones
CONTRIBUTING.md       the rules the oracle imposes; read before changing src/
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

## Licence

MIT; see `LICENSE`.

`reference/csegment/` is not MIT. It is the original IPW-based C, redistributed
under the University of California, Santa Barbara BSD licence, and it is not
needed to build or use the Rust program. `NOTICE` has the full terms. Its licence
asks that distributions including binaries carry this acknowledgement, so:

> This product includes software developed by the Computer Systems Laboratory,
> University of California, Santa Barbara and its contributors.

## Credits and citation

Original algorithm and C implementation by Jud Harward and Curtis Woodcock,
Boston University, building on the IPW library from UC Santa Barbara.

> Woodcock, C. and V. J. Harward. 1992. *Nested-hierarchical scene models and
> image segmentation.* International Journal of Remote Sensing 13(16): 3167-3187.

The second phase, and the parameters used here, are from:

> Ye, E., N. C. Coops, M. A. Wulder and T. Hermosilla. 2025. *A multi-resolution
> forest stand segmentation algorithm integrating Landsat imagery and forest
> structural, age, and species attributes.* ISPRS Journal of Photogrammetry and
> Remote Sensing. https://doi.org/10.1016/j.isprsjprs.2025.05.023

If you use this program in published work, please cite the Ye et al. paper for
the two-phase method and Woodcock & Harward for the algorithm underneath it.

TODO: add a LICENSE, and mint a DOI so the software itself is citable.
