# Temporal fixtures — segmentation over a time series

`tests/temporal/` is the oracle for `standseg temporal`, the multi-temporal
method of Ye, Coops, Hermosilla, Wulder & Gergel, "Temporally Consistent Forest Stand Segmentation Using Landsat Imagery", *IEEE Geosci. Remote Sens. Lett.* 22 (2025) 2506005, [doi:10.1109/LGRS.2025.3602095](https://doi.org/10.1109/LGRS.2025.3602095):
segment the first image of a series with both stages, then at every later step
re-segment only the pixels whose class changed, paste those segments into the
previous result, and develop the combined map against that step's second-stage
layer.

| | `tests/golden/` | `tests/stage2/` | `tests/temporal/` |
|---|---|---|---|
| Producer | the 1992 C | Elaine's Python stage 2 | Elaine's paper-2 pipeline (`forest-stand-update`), re-run here |
| Status | **Inviolable** | Regenerable, deliberately | Regenerable, deliberately |
| Verify | `verify_golden.sh` | `verify_stage2.sh` | `verify_temporal.sh` |

Regenerate only with `tools/temporal_oracle/gen_fixtures.py`, which writes to
`build/out/temporal/fixtures`; review that, then unlock, copy, re-pin, re-lock.
`tools/` is gitignored because the upstream repository is private.

## The oracle

`tools/temporal_oracle/pipeline.py` drives `run_our_algorithm.py` step for step
with the upstream helpers imported **verbatim** (`forest-stand-update` @
`9912421`; the algorithm last changed at `037a414`). Only the two external
programs it shells out to are replaced:

- **Stage 1** was the 1992 C on a Linux box; here it is `standseg`, which
  reproduces it byte for byte. The C never looked at nodata, but the
  normalized inputs declare `data ignore value = 0`, so the oracle strips that
  line before stage 1 sees it. With it honoured, FIRE1-1984 takes 29 passes
  instead of the archive's 31.
- **Stage 2** was `segment_python.stage2` @ `009a1ff`; here it is the vendored
  `tools/stage2_oracle` with `KEEP_ZERO=True` (the one behavioural difference
  between the two commits, shown to change no map) and the canonical tie rule.

Two guards abort the run, and neither has fired on any pinned data or on the
three full sites. One is a combined map with no 0 pixel, where upstream would
renumber a real segment to 0. The other is an id count past the inherited
dtype, where upstream would wrap.

How far the oracle is checked against the archived 2024 runs on the drive
(`validate_archive.py`, FIRE1/2/3 × 39 years) is shown below. The archive
predates the combine fix (`981892d`), so the old combine is checked verbatim.

| piece | result |
|---|---|
| normalized stage-1 and stage-2 inputs | byte-identical, 117/117 each |
| change masks | byte-identical, 114/114 |
| masked stage-1 maps, rmap + armap + pass counts | byte-identical, 117/117 |
| pre-fix combine on the archive's own maps | byte-identical, 114/114 |
| stage 2 | **not** reproducible; see PLAN.md §14.4 |

## The four cases

| case | grid | steps | source | tie sweep (distinct maps per step) |
|---|---|---|---|---|
| `tiny_synthetic` | 12×12×1 | 1–3 | hand-made | 1, 1, 1 |
| `fire1_crop` | 80×80×6 | 2001–2007 | FIRE1 (11S) rows 344–423, cols 224–303 | 2, 1, 1, 1, 1, 1, 1 |
| `fire2_crop` | 80×80×6 | 1996–2002 | FIRE2 (10S) rows 152–231, cols 192–271 | 4, 4, 4, 2, 2, 3, 3 |
| `fire3_crop` | 80×80×6 | 2001–2007 | FIRE3 (10S) rows 112–191, cols 712–791 | 1 × 7 |

The crops were chosen for change at every step and for some layer nodata, so
that the combine has its 0. The windows sit over each site's burn:

- fire1: 2,925 px change in 2003.
- fire2: 2,853 px change in 1999.
- fire3: 6,215 of 6,400 px change in 2004.

The three crops use the paper's parameters: `-t 50 -m 0.2 -n 18,18,18,36,36`,
stage 2 `80,8000`. `tiny_synthetic` uses `-t 50 -m 1 -n 16,16,32,200,200`,
stage 2 `6,200`. It has a left half, a right half and a nodata corner. Step 2
pastes a 4×4 square into the left half. Step 3 pastes a stripe at column 3
that cuts a 4-px piece off the left segment and another off the square. The
combine must number those pieces separately (raster order of each piece's
first pixel within its source id), and stage 2 must then fold both 4-px pieces
into their nearest neighbour in the layer. Every value in it can be checked by
hand.

The tie sweep runs every case under all six (candidate order × near-tie)
variants. It only records which steps are tie-sensitive; the pinned map is
always the canonical rule (ascending id, keep the incumbent), as in
`tests/stage2`.

### Layout

```
<case>/case.json                  parameters, source window, per-step counters and sweep
<case>/input/manifest.csv         label,image,stage2,classes -- what the CLI reads
<case>/input/<label>/image        raw stage-1 image (uint16, nodata 65535)
<case>/input/<label>/stage2       raw second-stage layer (uint16, nodata 65535)
<case>/input/<label>/classes      class layer (uint8, nodata 255)
<case>/expected/<label>.rmap.<n>  stage-1 map, n = stage-1 passes
<case>/expected/<label>.change    change mask (steps after the first)
<case>/expected/<label>.combined  combined map (steps after the first)
<case>/expected/<label>.armap.<n> final map, n = stage-2 passes
```

Inputs are the **raw** files, so the fixtures pin the paper's normalization
(`--normalize`) as well. All are ENVI, BSQ, little-endian, with a minimal
`.hdr`. The whole directory is marked binary in `.gitattributes`, so no
checkout can rewrite a byte.

`tests/temporal_match.rs` runs every case through the library and through the
CLI (from the manifest), comparing all 88 expected files byte for byte with the
pass counts. It also covers:

- `--change` masks against `--classes`;
- a negative control;
- a step with no change.
