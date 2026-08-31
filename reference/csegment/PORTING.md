# The C reference — what was changed to build it on macOS

This is the original 1990s C, vendored from `~/mac2025/segment/segment`, adapted
only as far as necessary to build and run on macOS/arm64 with clang. It exists as
a **debugging oracle**: something to instrument and perturb when the Rust
diverges. It is not the deliverable.

Verified: reproduces `proof/regmap.rmap.51` and `proof/regmap.armap.58`
byte-exactly, with `myseg.log` matching on every numeric value across all
51 normal + 58 auxiliary passes.

    cd reference/csegment && make
    cd reference/out/verify
    cp ../../../tests/golden/misc/temp_byte_bip{,.hdr} . && chmod 644 temp_byte_bip*
    ../../csegment/bin/segment -t 10 -m .1 -n 15,15,100,2500,2500 -o t10 temp_byte_bip

## Algorithm files

- `src/pixel.c`, `src/segment.c` — **untouched**.
- `src/region.c` — one change: `flip()` is wrapped in `#ifdef FORCE_FLIP` so the
  tie-break experiment can force it to a constant. A no-op unless `-DFORCE_FLIP`
  is passed.
- `src/set.c` — **the fix that actually mattered.** `add_to_set` read a 4-byte
  `REGION_ID` through a `long *` and deduped on all 8 bytes (the original source
  flags this itself: `// OFFENDER`). On Linux/x86-64 the adjacent stack garbage
  was stable within a `reg_nnbr` call, so it behaved as a 32-bit compare. On
  macOS/arm64 it is not, duplicate ids escape dedup, and every duplicate creates
  an exact-distance tie that burns an extra `flip()` draw. `case 4` is now
  well-defined 32-bit; `-DSET_LEGACY_LONG_UB` restores the literal original,
  which diverges.

## Support files

- `src/gdal_io.c` — replaced (original kept as `.gdal-orig`). GDAL is not
  installed and installing it was out of scope; this reads raw ENVI uint8 BIP via
  the `.hdr` sidecar and writes the region map with the same `nbits` ladder and
  byte layout the GDAL ENVI driver produced. IPW input is **not** supported here,
  so Case 2 cannot be run against this oracle.
- `src/glibc_random.c` — a glibc-compatible `random()`. Kept as a guarantee
  against libc drift, but **not required**: macOS `random()` turns out to emit the
  identical sequence, and the build reproduces the golden without it.
- `src/linux_compat.c` — defines `getpagesize()` as 4096. This is load-bearing
  for a line-for-line log match: `reclaim_trigger` derives from it (911 on Linux,
  3641 on Apple silicon), which changes *when* the region list is compacted.
  Payloads are byte-exact either way, so compaction is output-neutral.
- `inc/values.h` — shim for `MAXSHORT`/`MAXINT`/`MAXLONG`/`MAXFLOAT`.
- `inc/extern.h` — uses the system `<string.h>` instead of a vendored glibc
  header; `sys_nerr` disabled (macOS declares it `const`, nothing uses it).
- `Makefile` — original kept as `.orig`. No `gdal-config`; `-O0
  -ffp-contract=off` plus warning suppressions for K&R-era constructs.

## Licence cleanup

Two files under `inc/PORT/linux/` were bundled copies of GCC headers, GPLv2 with
the header exception, and the only GPL-licensed code in a repository that is
otherwise MIT and BSD. Both are gone.

- `string.h` was never referenced. No `ALT_STRING_H` is defined anywhere.
- `float.h` was reached through `#define ALT_FLOAT_H` in `PORT/linux/config.h`,
  which `inc/config.h` symlinks to. `ipw.h` says why it exists: *"if
  /usr/include/float.h is missing, include local float.h in config.h"* — a
  fallback for pre-standard systems. With `ALT_FLOAT_H` left undefined, `ipw.h`
  falls through to the system `<float.h>`.

Checked, not assumed: after a clean rebuild the C still reproduces both Case 1
golden maps byte for byte (`regmap.rmap.51` and `regmap.armap.58`).

## Provenance note

Commits `0c1e9c4` through `b584723` describe Rust work but also carry snapshots
of this C porting effort, which was running concurrently in the same working
tree and was swept in by `git add -A`. The C changes are documented here rather
than in those commit messages.
