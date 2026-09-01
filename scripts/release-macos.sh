#!/usr/bin/env bash
#
# Build -> smoke-test -> sign -> notarize -> package the macOS binary.
#
# macOS is not in .github/workflows/release.yml because signing needs a
# Developer ID certificate, and putting one in GitHub Secrets is a bigger
# commitment than this project needs. So the macOS artifact is built here, on a
# Mac with the certificate already in its keychain, and uploaded by hand.
#
# Output: build/out/release/standseg-universal-apple-darwin.tar.gz
#
# Prereqs (one-time):
#   - A "Developer ID Application" certificate in the login keychain.
#   - A notarytool keychain profile holding an Apple ID + app-specific password:
#       xcrun notarytool store-credentials standseg-notary \
#           --apple-id <you@example.com> --team-id <TEAMID> --password <app-specific>
#     The profile is account credentials, not per-app, so an existing profile
#     from another project works -- pass it as NOTARY_PROFILE.
#
# Usage:
#   ./scripts/release-macos.sh
#   NOTARY_PROFILE=other-profile ./scripts/release-macos.sh
#   ./scripts/release-macos.sh v0.2.0     # also uploads to that GitHub release
#
# Nothing here is credential-bearing: the identity is looked up in the keychain
# at run time and the notary password never leaves it.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO/build/out/release"
NOTARY_PROFILE="${NOTARY_PROFILE:-standseg-notary}"
TAG="${1:-}"
ARCHS=(aarch64-apple-darwin x86_64-apple-darwin)

step() { printf "\n\033[1;34m==> %s\033[0m\n" "$*"; }
ok()   { printf "\033[1;32m  ok  %s\033[0m\n" "$*"; }
die()  { printf "\033[1;31mERROR: %s\033[0m\n" "$*" >&2; exit 1; }

cd "$REPO"

step "Locating the signing identity"
IDENTITY="${SIGN_IDENTITY:-$(security find-identity -v -p codesigning \
    | sed -n 's/.*"\(Developer ID Application: .*\)"/\1/p' | head -n1)}"
[[ -n "$IDENTITY" ]] || die "no Developer ID Application certificate in the keychain"
# Print only the type, never the team id -- this runs in a terminal that gets
# pasted into issues.
ok "Developer ID Application certificate found"

xcrun notarytool history --keychain-profile "$NOTARY_PROFILE" >/dev/null 2>&1 \
    || die "no notarytool keychain profile '$NOTARY_PROFILE' (see the header of this script)"
ok "Notary profile '$NOTARY_PROFILE' works"

step "Building ${#ARCHS[@]} architectures"
for a in "${ARCHS[@]}"; do rustup target add "$a" >/dev/null; done
for a in "${ARCHS[@]}"; do
  # Stripped, matching the Linux and Windows artifacts. `debug = 1` stays in the
  # profile for local work; a download does not need the symbols.
  CARGO_PROFILE_RELEASE_STRIP=symbols cargo build --release --locked --target "$a"
  ok "$a"
done

step "Merging into one universal binary"
rm -rf "$OUT" && mkdir -p "$OUT"
BIN="$OUT/standseg"
lipo -create -output "$BIN" \
    "target/aarch64-apple-darwin/release/standseg" \
    "target/x86_64-apple-darwin/release/standseg"
lipo -info "$BIN"

# Same check the CI release job runs: segment the reference scene with the
# binary that is about to ship and require the 1992 C's exact bytes. Output goes
# to build/out, never into tests/.
step "Segmenting the reference scene with the built binary"
"$BIN" --version
"$BIN" -t 10 -m .1 -n 15,15,100,2500,2500 \
    -o demo --outdir "$OUT/smoke" tests/golden/misc/temp_byte_bip >/dev/null
cmp "$OUT/smoke/demo.rmap.51"  tests/golden/test_3456/expected/proof/regmap.rmap.51
cmp "$OUT/smoke/demo.armap.58" tests/golden/test_3456/expected/proof/regmap.armap.58
rm -rf "$OUT/smoke"
ok "region maps are byte-identical to the reference output"

step "Signing"
# --options runtime is the hardened runtime, which notarization requires.
# --timestamp gets a trusted timestamp so the signature outlives the cert.
codesign --force --options runtime --timestamp --sign "$IDENTITY" "$BIN"
codesign --verify --strict --verbose=2 "$BIN"
ok "Signed with a hardened runtime"

step "Notarizing (2-10 min)"
ZIP="$OUT/notarize.zip"
ditto -c -k --keepParent "$BIN" "$ZIP"
xcrun notarytool submit "$ZIP" --keychain-profile "$NOTARY_PROFILE" --wait
rm -f "$ZIP"
ok "Notarized"

# No `stapler staple` here, and that is not an oversight: a staple can only be
# attached to a bundle, a .dmg or a .pkg, never to a bare Mach-O executable.
# Gatekeeper therefore confirms this binary online the first time it runs, which
# needs a network connection once. Shipping a .dmg purely to allow stapling
# would be a strange thing to hand someone for a command line tool.
#
# `spctl --assess --type exec` reports "rejected (the code is valid but does not
# seem to be an app)" for a signed, notarized CLI binary. That is spctl refusing
# to judge a non-bundle, not Gatekeeper refusing to run it, so it is reported
# here and not treated as a failure. The check that decides whether a download
# works is at the end of this script.
step "Gatekeeper assessment (informational)"
spctl --assess --type exec -vv "$BIN" 2>&1 | sed -E 's/Developer ID Application: [^(]*/Developer ID Application: /' || true

step "Packaging"
D="$OUT/standseg-universal-apple-darwin"
mkdir -p "$D"
mv "$BIN" "$D/"
cp LICENSE README.md "$D/"
tar czf "$D.tar.gz" -C "$OUT" "$(basename "$D")"
rm -rf "$D"
ok "$D.tar.gz"

# The real test: unpack the artifact somewhere else, stamp it with the
# quarantine attribute a browser applies to a download, and run it. If signing
# or notarization had gone wrong this is where macOS refuses, with "cannot be
# opened because the developer cannot be verified".
step "Running the packaged binary as a downloader would"
Q="$(mktemp -d)"
tar xzf "$D.tar.gz" -C "$Q"
QBIN="$Q/$(basename "$D")/standseg"
xattr -w com.apple.quarantine "0083;$(printf %x "$(date +%s)");Safari;$(uuidgen)" "$QBIN"
"$QBIN" --version >/dev/null || die "a quarantined copy will not run -- do not ship this"
rm -rf "$Q"
ok "a quarantined copy runs; Gatekeeper is satisfied"

if [[ -n "$TAG" ]]; then
  step "Uploading to release $TAG"
  gh release upload "$TAG" "$D.tar.gz" --clobber
  ok "Uploaded"
else
  step "Done"
  echo "  Upload it with:  gh release upload <tag> '$D.tar.gz'"
fi
