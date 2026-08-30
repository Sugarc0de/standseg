#!/usr/bin/env bash
# Verify the stage-2 fixtures are byte-identical to what tools/stage2_oracle
# produced. Exits non-zero (and names the offenders) on any drift.
#
# Unlike tests/golden, these ARE regenerable -- but only deliberately, by running
# the generator and reviewing the diff. Drift detected here means something wrote
# into tests/stage2, which is exactly what must not happen.
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$here/stage2" || { echo "tests/stage2 missing" >&2; exit 1; }

actual=$(find . -type f | sed 's|^\./||' | sort)
expect=$(sed 's/^[0-9a-f]*  //' "$here/STAGE2.sha256" | sort)
if [ "$actual" != "$expect" ]; then
    echo "STAGE-2 FIXTURE SET CHANGED (files added or removed):" >&2
    diff <(echo "$expect") <(echo "$actual") | sed 's/^/  /' >&2
    exit 1
fi

if ! shasum -a 256 -c "$here/STAGE2.sha256" --quiet; then
    echo "STAGE-2 FIXTURES MODIFIED -- restore them from git before continuing." >&2
    echo "Regenerate only via tools/stage2_oracle/gen_fixtures.py, never by hand." >&2
    exit 1
fi

echo "stage-2 fixtures OK ($(wc -l < "$here/STAGE2.sha256" | tr -d ' ') files)"
