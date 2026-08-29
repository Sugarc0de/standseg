# Golden fixtures — provenance and how to compare

Captured from `~/mac2025/segment/dataset` (Harward & Woodcock segmentation, IPW-era
reference outputs). `tests/golden/` is **read-only and checksum-pinned**; it is the
only oracle for the rewrite. Never edit, regenerate, or "fix" anything under it.

Verify at any time: `tests/verify_golden.sh`

## Invocation recorded in the fixtures

Every expected output carries the original command in its IPW header:

    segment -t 10 -m .1 -n 15,15,100,2500,2500 -o t10-m1-n15_15_100_2500_2500_myseg <input>.ipw

Same parameters for both cases (this is also what `misc/param.txt` holds).

## Case 1 — `test_3456` (primary)

| | |
|---|---|
| Input | `test_3456/input/test_3456.bip` (ENVI) / `test_3456.bip.ipw` (IPW) |
| Geometry | 250 × 250, 4 bands, ENVI data type 1 (uint8), BIP |
| Converged | region map 51 passes, aux region map 58 passes |

Expected outputs in `test_3456/expected/`:

- `..._myseg.armap.58` — final aux region map, **IPW container**: 320-byte text header + 125000 bytes payload
- `..._myseg.rmap.51` — final region map, IPW container, same layout
- `regionmap`, `regionmap.hdr`, `regionmap.aux.xml` — ENVI copy of the final armap
- `proof/regmap.armap.58`, `proof/regmap.rmap.51` — raw ENVI payloads (125000 bytes, no header) + `.hdr` sidecars
- `myseg.log` — full per-pass merge statistics

**The 320-byte size difference is a container difference, not a data difference.**
Verified: the IPW payload is byte-identical to the `proof/` files.

    tail -c 125000 expected/..._myseg.armap.58  ==  proof/regmap.armap.58   (sha256 3ed2dbdc…)
    tail -c 125000 expected/..._myseg.rmap.51   ==  proof/regmap.rmap.51    (sha256 1d79a250…)
    expected/regionmap                          ==  proof/regmap.armap.58

So `proof/` and the IPW files are the same oracle in two wrappers. Compare against
whichever matches the container your rewrite emits.

## Case 2 — `LC80220492014083LGN00` (secondary)

| | |
|---|---|
| Input | `input/LC80220492014083LGN00_stack` (ENVI) / `_stack.ipw` (IPW) |
| Geometry | 250 × 250, 8 bands, ENVI data type 2 (int16), BIP |
| Converged | region map 17 passes, aux region map 1 pass |

Aux segmentation completed in a single pass here, so `armap.1`, `rmap.17` (payload)
and `regionmap` are all byte-identical (sha256 `fad51722…`). IPW header is 333 bytes.

## `misc/`

`test.bip` (62500 B), `temp_byte_bip` (250000 B) and their sidecars. No expected
output was committed for these — usable as smoke inputs, not as an oracle.

## Comparing a rewrite run

Write program output to `build/out/` (gitignored). Never point the program at
`tests/golden/` — the original repo kept inputs and outputs in one directory, which
is exactly how the oracle gets clobbered.

    # payload-only comparison, container-agnostic
    cmp <(tail -c 125000 build/out/..._myseg.armap.58) tests/golden/test_3456/expected/proof/regmap.armap.58

    # first differing byte offset, if any
    cmp -l <(...) <(...) | head

Byte-exact equality on the payload is the pass condition. `myseg.log` is useful for
localising a divergence (which pass, which merge count) but is not itself a pass/fail
criterion — it contains timings and paths.
