# fast_segment

Read `CONTRIBUTING.md` first and follow it. It carries the rules that matter
here: the golden fixtures are inviolable and must never be edited, regenerated,
or re-baselined; nothing is ever written into `tests/`; all output goes to
`build/out/`; passing means a byte-exact region-map payload, never "close
enough"; and commit messages state the real comparison result.

Two things specific to working here with an agent:

- If a fixture looks wrong, stop and ask rather than working around it. A
  mismatch is a bug in this program until proven otherwise.
- Never run `tests/lock_golden.sh --unlock`, and never push.

`PLAN.md` has the design notes, the port hazards (§3), and the second phase
(§13). The original C this reimplements lives at `~/mac2025/segment`; it is
deliberately not on the tool path.
