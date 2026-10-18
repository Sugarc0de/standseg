#!/usr/bin/env bash
# Verify the temporal fixtures are byte-identical to what tools/temporal_oracle
# produced. Exits non-zero (and names the offenders) on any drift.
#
# Unlike tests/golden, these ARE regenerable -- but only deliberately, by running
# the generator and reviewing the diff. Drift detected here means something wrote
# into tests/temporal, which is exactly what must not happen.
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$here/temporal" || { echo "tests/temporal missing" >&2; exit 1; }

actual=$(find . -type f | sed 's|^\./||' | sort)
expect=$(sed 's/^[0-9a-f]*  //' "$here/TEMPORAL.sha256" | sort)
if [ "$actual" != "$expect" ]; then
    echo "TEMPORAL FIXTURE SET CHANGED (files added or removed):" >&2
    diff <(echo "$expect") <(echo "$actual") | sed 's/^/  /' >&2
    exit 1
fi

if ! shasum -a 256 -c "$here/TEMPORAL.sha256" --quiet; then
    echo "TEMPORAL FIXTURES MODIFIED -- restore them from git before continuing." >&2
    echo "Regenerate only via tools/temporal_oracle/gen_fixtures.py, never by hand." >&2
    exit 1
fi

echo "temporal fixtures OK ($(wc -l < "$here/TEMPORAL.sha256" | tr -d ' ') files)"
