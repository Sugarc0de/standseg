# fast_segment — rewrite rules

A from-scratch reimplementation of the Harward & Woodcock nested-hierarchical image
segmentation algorithm (see `~/mac2025/segment`, not on the tool path — do not add it).
The reference outputs it produced are vendored here as the oracle.

## The oracle is inviolable

`tests/golden/` and `tests/GOLDEN.sha256` are the captured reference outputs. They are
chmod 444/555 and checksum-pinned, and a `Stop` hook verifies them after every turn.

- Never edit, delete, move, or regenerate anything under `tests/golden/`.
- Never regenerate `tests/GOLDEN.sha256`. If verification fails, the fixtures are
  damaged — restore them from git, don't re-baseline them.
- Never run `tests/lock_golden.sh --unlock`.
- If a fixture looks wrong, **stop and ask**. A mismatch is a bug in the rewrite until
  proven otherwise.

## Never write into tests/

The original repo kept program inputs and outputs in the same directory, which is how
an oracle gets clobbered by its own program. Here:

- All program output goes to `build/out/` (gitignored).
- Inputs are read from `tests/golden/*/input/` and never written back.
- Read `tests/GOLDEN.md` before writing any comparison code — it documents the two
  cases, the IPW-vs-ENVI container difference, and the exact comparison recipe.

## Pass condition

Byte-exact equality of the region-map payload against the golden output, for the
invocation recorded in the fixture headers:

    -t 10 -m .1 -n 15,15,100,2500,2500

`myseg.log` is a debugging aid for localising divergence, not a pass/fail criterion.
"Close enough", "visually similar", or "differs in N pixels" is a failure. Report the
actual byte comparison, never an assertion that it passed.

## Commits

- Commit after each unit of work that builds and leaves the golden check passing.
- Never push. Never amend, reset, revert, or rebase — the history is the safety net.
- Commit messages state what was implemented and the real comparison result, e.g.
  `armap payload matches golden (125000/125000 bytes)` or `still diverges at pass 12`.
  Never claim a match without having run the comparison in that same session.
