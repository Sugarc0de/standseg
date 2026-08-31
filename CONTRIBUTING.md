# Contributing

The unusual thing about this repository is that it has an oracle. Almost every
change can be checked against a program that already produced the right answer,
and the project is built around not giving that up. Most of the rules below
follow from it.

## The golden fixtures are inviolable

`tests/golden/` holds the reference inputs and outputs captured from the original
C. `tests/stage2/` holds cases generated from the reference Python. Both are
checksum-pinned (`tests/GOLDEN.sha256`, `tests/STAGE2.sha256`) and stored
read-only.

- Never edit, delete, move, or regenerate anything under `tests/golden/`.
- Never regenerate `tests/GOLDEN.sha256`. If verification fails, the fixtures are
  damaged. Restore them from git; do not re-baseline them.
- If a fixture looks wrong, stop and ask. A mismatch is a bug in this program
  until proven otherwise, and that has been true every time so far.

Git does not preserve read-only permissions, so after cloning:

```bash
tests/lock_golden.sh      # make the fixtures read-only
tests/verify_golden.sh    # confirm nothing has drifted
```

`tests/stage2/` differs in one way: it *is* regenerable, from the stage-2
Python. `tests/golden/` never is.

## Never write into tests/

The original program kept inputs and outputs in the same directory, which is how
an oracle gets clobbered by the program it is meant to check.

- All program output goes to `build/out/`, which is gitignored.
- Inputs are read from `tests/golden/*/input/` and never written back.
- Read `tests/GOLDEN.md` before writing comparison code. It documents the two
  cases, the IPW-versus-ENVI container difference, and the exact recipe.

## What counts as passing

Byte-exact equality of the region-map payload against the reference output, for
the invocation recorded in the fixture headers:

    -t 10 -m .1 -n 15,15,100,2500,2500

"Close enough", "visually similar", and "differs in N pixels" are failures. Two
segmentations can look identical and disagree everywhere that matters, so the
only useful report is the byte comparison itself. Quote the actual result rather
than asserting that it passed.

`myseg.log` is a debugging aid for localising a divergence, not a pass criterion.

## Running the checks

```bash
cargo test --release                        # everything, including both oracles
cd tests/golden && shasum -a 256 -c ../GOLDEN.sha256
bash tests/verify_stage2.sh
```

CI runs the same on Linux, macOS and Windows.

`cargo fmt` and `cargo clippy` are advisory. The code predates any formatter and
contains a few deliberate offences -- `out[0 * 9 + 0]` is `row * width + col`,
which reads better than the constant it folds to. Cleaning that up is welcome as
its own change; please do not mix it with a behavioural one, because a
reformatting diff hides exactly the kind of detail this project depends on.

## Commits

- Commit after each unit of work that builds and leaves the golden check passing.
- Commit messages should state what was implemented and the real comparison
  result, e.g. `armap payload matches golden (125000/125000 bytes)` or `still
  diverges at pass 12`. Never claim a match without having run the comparison.
- Do not amend, reset, revert, or rebase. The history is the safety net.

## Where things are

```
src/                  the segmenter
tests/golden/         1992 reference inputs and outputs, checksum-pinned
tests/stage2/         two-image segment-development fixtures, checksum-pinned
PLAN.md               design notes: algorithm, port hazards, memory, milestones
```

The two oracles are **not in this repository**. The 1992 C lives at
`~/mac2025/segment` and the stage-2 Python at `~/mac2025/segment_python`
(`tests/STAGE2.md` names the commit). They are someone else's code and thesis
code respectively, and neither is needed to build, run or test this program --
but they are still what a byte comparison is run against, so anything below that
says "against the C" or "against the Python" means those working copies. Docs
that name `reference/csegment/` or `tools/stage2_oracle/` mean a local checkout
at those paths; both are gitignored.

`PLAN.md` is worth reading before changing anything in `src/`. Section 3 lists
the details of the original that are easy to get wrong, and section 13 covers the
second phase.
