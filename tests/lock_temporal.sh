#!/usr/bin/env bash
# Make tests/temporal read-only (default) or writable (--unlock).
# Unlock only to regenerate fixtures from tools/temporal_oracle; re-lock after.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" = "--unlock" ]; then
    chmod -R u+w "$here/temporal"
    chmod u+w "$here/TEMPORAL.sha256"
    echo "tests/temporal UNLOCKED -- re-lock with: tests/lock_temporal.sh"
else
    find "$here/temporal" -type f -exec chmod 444 {} +
    find "$here/temporal" -type d -exec chmod 555 {} +
    chmod 444 "$here/TEMPORAL.sha256"
    echo "tests/temporal LOCKED (read-only)"
fi
