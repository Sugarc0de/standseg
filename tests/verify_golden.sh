#!/usr/bin/env bash
# Verify the golden fixtures are byte-identical to what was captured from the
# original `segment` repo. Exits non-zero (and names the offenders) on any drift.
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$here/golden" || { echo "tests/golden missing" >&2; exit 1; }

# Catch files that were added or removed as well as ones that changed.
actual=$(find . -type f | sed 's|^\./||' | sort)
expect=$(sed 's/^[0-9a-f]*  //' "$here/GOLDEN.sha256" | sort)
if [ "$actual" != "$expect" ]; then
    echo "GOLDEN FIXTURE SET CHANGED (files added or removed):" >&2
    diff <(echo "$expect") <(echo "$actual") | sed 's/^/  /' >&2
    exit 1
fi

if ! shasum -a 256 -c "$here/GOLDEN.sha256" --quiet; then
    echo "GOLDEN FIXTURES MODIFIED -- restore them before continuing." >&2
    echo "They are the only oracle for the rewrite; do not regenerate the manifest." >&2
    exit 1
fi

echo "golden fixtures OK ($(wc -l < "$here/GOLDEN.sha256" | tr -d ' ') files)"
