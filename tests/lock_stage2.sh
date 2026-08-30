#!/usr/bin/env bash
# Make tests/stage2 read-only (default) or writable (--unlock).
# Unlock only to regenerate fixtures from tools/stage2_oracle; re-lock after.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--unlock" ]; then
    chmod -R u+w "$here/stage2"
    chmod u+w "$here/STAGE2.sha256"
    echo "tests/stage2 UNLOCKED -- re-lock with: tests/lock_stage2.sh"
else
    find "$here/stage2" -type f -exec chmod 444 {} +
    find "$here/stage2" -type d -exec chmod 555 {} +
    chmod 444 "$here/STAGE2.sha256"
    echo "tests/stage2 LOCKED (read-only)"
fi
