#!/bin/bash
# Assembles the Phy CE web page in _site/: the page (site/), its screens (docs/phyce/, drawn with
# the app's own steps page: tools/hostgiac -V) and the bundle for the web installer.
#   site/build.sh <PhyCE.b84> [version]
# The workflow (.github/workflows/phyce.yml) builds the app, bundles it and runs this.
set -euo pipefail
cd "$(dirname "$0")/.."
B84=${1:?usage: site/build.sh <PhyCE.b84> [version]}
VER=${2:-}
OUT=_site
rm -rf "$OUT"
mkdir -p "$OUT/media/screenshots" "$OUT/files"
cp site/index.html site/style.css site/installer.js site/dusb.js "$OUT/"
cp -R site/artifice-v3 "$OUT/"
python3 site/fetch_inequalz.py "$OUT/artifice-v3" ${FLOWCE_INEQUALZ_FILE:+"$FLOWCE_INEQUALZ_FILE"}
cp site/icon.svg "$OUT/media/icon.svg"
cp docs/phyce/*.png "$OUT/media/screenshots/"
cp "$B84" "$OUT/files/PhyCE.b84"
printf '%s\n' "$VER" > "$OUT/files/version.txt"
echo "built $OUT ($(du -sh "$OUT" | cut -f1)), Phy CE ${VER:-(no version)}"
