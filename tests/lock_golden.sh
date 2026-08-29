#!/usr/bin/env bash
# Make tests/golden read-only (default) or writable (--unlock).
# Unlock only to intentionally re-capture fixtures; re-lock immediately after.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--unlock" ]; then
    chmod -R u+w "$here/golden"
    chmod u+w "$here/GOLDEN.sha256"
    echo "tests/golden UNLOCKED -- re-lock with: tests/lock_golden.sh"
else
    find "$here/golden" -type f -exec chmod 444 {} +
    find "$here/golden" -type d -exec chmod 555 {} +
    chmod 444 "$here/GOLDEN.sha256"
    echo "tests/golden LOCKED (read-only)"
fi
