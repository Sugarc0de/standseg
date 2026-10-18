# Checklist — multi-temporal segmentation (paper 2)

Live progress for adding the paper-2 algorithm (Ye et al., IEEE, doc 11137370):
segmentation **over a time series of images**, where each step re-segments only
what changed since the previous step and carries the rest forward. Same
discipline as paper 1 (`PLAN.md` §13): the reference Python is re-run here to
make an oracle, and the Rust has to reproduce it **byte for byte**, every time step.

Legend: `[x]` done and verified (the evidence is on the line) · `[~]` in
progress · `[ ]` not started · `[?]` needs Elaine.

---

## 0. What the algorithm is (recovered from `forest-stand-update`)

Source: `Sugarc0de/forest-stand-update`, `scripts/run_our_algorithm.py` +
`scripts/helpers/{image_handler,io_handler}.py`. The algorithm last changed at
`037a414` (2025-02-09); HEAD `9912421` touches only plots. Stage 2 is
`segment_python.stage2` at `009a1ff`.

For a sequence of time steps `t = 0 … T-1`:

| step | t = 0 (reference) | t ≥ 1 |
|---|---|---|
| prepare | rescale stage-1 and stage-2 images to 1..255 per band, nodata → 0 | same |
| change mask | — | `M_t = class_t != class_{t-1}` (1 = changed) |
| stage 1 | Woodcock-Harward normal passes on image_t → `U_0` | the same, but **only inside `M_t`** (`-M`) → `U_t` |
| combine | — | `C_t` = previous result with `U_t` pasted over the changed pixels, then every 4-connected piece of every id gets its own id, numbered 1..n |
| stage 2 | segment development of `U_0` on layer_0 → `R_0` | segment development of `C_t` on layer_t → `R_t` |

Paper parameters: `-t 50 -m 0.2 -n 18,18,18,36,36`, stage 2 `80,8000`, change
from VLCE 2.0 land cover, stage-2 layer = gross stem volume.

Nothing in it is forest-specific. "Class layer" is any categorical map whose
change between steps defines what gets re-segmented: land cover, a burn
severity class, a disturbance product. The Rust names follow that.

## 1. Study data on the drive

`/Volumes/easystore/UBC/second_project/paper_data/study_sites/` (snapshot of
2024-10-09/10, also in `second_project.zip`):

| site | UTM | grid | biggest change year | what's there |
|---|---|---|---|---|
| FIRE1 | 11S | 585 × 484 | 2004 (12.3 % of px) | proxies, VLCE, GSV/p95/cv 1984–2022 + full archived run `results/1984…2022` |
| FIRE2 | 10S | 446 × 446 | 1999 (5.0 %) | same |
| FIRE3 | 10S | 1056 × 736 | 2004 (24.9 %) | same |
| HAR1–3 | — | — | — | polygon only, no rasters |
| `archive_site/site1` | 10S | | | older 39-year run (pre-Oct 2024 code) |
| `stikine_fire` | 8S/9S/mosaic | | | inputs, no results |

- [x] Located the paper-2 study data and the archived runs (above).
- [?] **Only three of the six paper sites are on this drive.** The published
      six are FIRE3/7/11/12/14/15 (Cariboo, Thompson-Nicola, Central Coast,
      Fraser-Fort George, Bulkley-Nechako, Peace River), run in 2025 from
      `/mnt/e/cfs/...`. This drive's FIRE1/2/3 are the October 2024 set under
      older names. Need the 2025 study_sites dir to cover all six.

## 2. Hazards found before writing any code

- [x] **The C ignored nodata; the Rust reads it from the header.** The
      normalized inputs declare `data ignore value = 0`. With it honoured,
      FIRE1-1984 stage 1 takes 29 passes; with it stripped, 31 — and
      `rmap.31`/`armap.19` are byte-identical to the 2024 Linux C output. The
      temporal command must not auto-apply header nodata to stage 1.
- [x] **Stage 1 under the change mask already reproduces the 2024 C**:
      FIRE1-1985 `-M` run, `rmap.11` and `armap.4` byte-identical.
- [x] **Python stage 2 changed after the vendored commit.** `009a1ff` (what
      paper 2 imported) builds region 0 like any other region; `427a5a3`
      (vendored in `tools/stage2_oracle`) skips it. Region 0 can never be a
      neighbour (adjacency excludes 0 on both sides), so maps should be
      identical and only the `considered`/`no_cand` counters differ.
      FIRE1-1984: identical (no region 0 there). To verify on years with a
      region 0.
- [x] **The archived runs predate the combine fix.** Archive is 2024-10-09;
      the connected-component relabel in `combine_two_region_map_with_no_change_mask`
      landed 2024-11-17 (`981892d`). So archived years ≥ 1985 used the old
      combine and are not the paper's final algorithm.
- [?] **The archived stage-2 maps are not reproducible from the committed
      Python.** FIRE1-1984: archive vs. re-run differ in the *zero set* by 389
      px — 346 regions of 1–2 px that are 100 % nodata were absorbed by a
      neighbour instead of dropped. The committed code cannot do that (it would
      drop them, and asserts on a zero centroid). Rust and the re-run Python
      agree exactly (0 px). Suggests the Linux box had uncommitted stage-2
      changes. Not blocking: the oracle is the committed code, as in paper 1.
- [ ] Unseeded tie-break in Python stage 2 → pin the canonical rule (ascending
      id, keep incumbent), as paper 1. Sweep the 6 variants per fixture to
      record which cases are tie-sensitive.
- [ ] `combine` renumbers with `enumerate(np.unique(updated_map))`: if the map
      has **no** zero pixel, the smallest real region becomes id 0 (silently
      masked). Guard in the oracle with an assert; Rust keeps 0 for 0 only.
      Expect it never to fire on real data — to confirm.
- [ ] `combine` writes with the previous map's dtype: uint16 could wrap past
      65 535 ids. Assert in oracle; Rust widens instead.
- [ ] Multiple `-t` tolerances: the Python takes `os.listdir(...)[0]` of
      several `.rmap.*` files — arbitrary. Rust uses the last (final) map.
      Paper used one tolerance, so not exercised.

## 3. Oracle (`tools/temporal_oracle/`)

- [x] `tools/temporal_oracle/pipeline.py` (local-only, like `tools/stage2_oracle`:
      `tools/` is gitignored and the upstream repo is private). Imports the
      upstream helpers **verbatim** from `upstream/` (@ 9912421); stage 1 =
      `standseg` with header nodata stripped; stage 2 = vendored Python with
      `KEEP_ZERO=True` (new switch in `tools/stage2_oracle/segment.py`, default
      off so `tests/stage2` is untouched), canonical tie rule. Guards abort on
      the no-zero and dtype-wrap hazards.
- [ ] Validate pieces against the archive, all 3 sites × 39 years:
  - [ ] normalization → `processed_input/*_norm.tif` byte-identical
  - [ ] change mask → `processed_input/*__mask.tif` byte-identical
  - [ ] stage 1 (masked) → archived `stage1.rmap.N` / `stage1.armap.N` byte-identical
  - [ ] old combine (pre-`981892d`) → archived `new_stage1_map.tif` byte-identical
  - [ ] stage 2 vs archive: quantify (zero-set + pixel residual), as §13.7
- [ ] `KEEP_ZERO` on vs off: identical maps on every year of every site.
- [~] Full-series reference outputs, 3 sites × 39 years → `build/out/temporal/oracle/`.
      FIRE1, FIRE2 done (39 steps each, guards never fired); FIRE3 running.
- [ ] Fixtures `tests/temporal/`: raw crops over each fire + one hand-checkable
      synthetic case, checksum-pinned (`TEMPORAL.sha256`), locked, with the
      tie-rule sweep recorded per case.

## 4. Rust

- [x] `src/temporal.rs` — generic, no forest wording:
  - [x] `normalize` (per-band rescale to 1..255, nodata → 0, the paper's prep)
  - [x] `change_mask` (class layer t vs t-1)
  - [x] `combine` (paste + 4-connected relabel, scipy numbering — checked
        empirically: scipy labels in raster order of first pixel, 3000/3000
        random masks)
  - [x] `Series::step`: step 0 = stage 1 + stage 2; step t = masked stage 1 → combine → stage 2
- [x] `driver::run_normal` (phase 0 + normal passes); `run_with_stage2` now
      calls it. `cargo test --release` all green, golden 31 files + stage-2 43
      files verify.
- [x] CLI `standseg temporal` (`--image`/`--stage2`/`--classes`|`--change`
      per step, or `--manifest`; `--normalize`; `--intermediates`). Dispatched
      before the flat parser, so existing invocations are untouched.
- [ ] Unit tests: combine numbering vs scipy, normalization vs archived files.
- [ ] Fixture tests (library + CLI), with a negative control.

## 5. Verification

- [ ] All fixtures byte-exact, every step.
- [~] Full series, FIRE1/2/3 × 1984–2022 = 117 steps, byte-exact against the oracle.
  - [x] FIRE1: 39/39 steps, every artifact (stage-1 map, change mask,
        combined map, final map) identical — 76 164 660 bytes, 0 differ.
  - [x] FIRE2: 39/39 steps — 53 707 320 bytes, 0 differ.
  - [ ] FIRE3
- [ ] Golden + stage-2 suites unchanged (`cargo test --release`, both manifests).
- [ ] Timing vs the Python pipeline.

## 6. Docs

- [ ] `PLAN.md` §14, `tests/TEMPORAL.md`, README usage section.
