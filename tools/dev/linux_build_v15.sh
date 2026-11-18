#!/bin/bash
# Phy CE: tools/dev/wsl_build_v15.sh for Linux (the workflow, or a Linux shell), building in place.
#   tools/dev/linux_build_v15.sh <repo> "APPLANG=en PHYS=1"   (CE toolchain v15 in ~/CEdev; APP_NAME, APP_INFO in the environment)
# Linux adaptation of tools/dev/wsl_build_v15.sh, building in place.
set -uo pipefail
SRC="$1"; ARGS="$2"
CEDEV15="$HOME/CEdev"; BU="$CEDEV15/binutils/bin"; export PATH="$CEDEV15/bin:$PATH"
cd "$SRC"
LDX="--defsym=___heapbot=___heap_low --defsym=___heaptop=___heap_high"
python3 tools/catpack.py src/menufr.cc src/catalog_text.h || exit 1
PAD=0
link() {
  rm -f bin/DEMO.obj
  make -j"$(nproc)" $ARGS APPLICATION=YES LINKER_SCRIPT="$CEDEV15/meta/linker_script_app.ld" EXTRA_LDFLAGS="$LDX" OUTPUT_MAP=YES STRIP_CMD=: LOAD_ADDR="$1" bin/DEMO.obj > build.log 2>&1 \
    || { echo "[build] FAILED"; grep -E "error|undefined" build.log | head -20; exit 1; }
  "$BU/z80-none-elf-objcopy" -O binary bin/DEMO.obj bin/body.bin
  [ "$PAD" = 0 ] || head -c "$PAD" /dev/zero | LC_ALL=C tr '\0' '\377' >> bin/body.bin
}
pack_all() {
  LA=$(python3 app_tools/pack_app.py size bin/body.bin); n1=$(stat -c %s bin/body.bin)
  link "$LA"; n2=$(stat -c %s bin/body.bin); [ "$n1" = "$n2" ] || { echo "size changed"; exit 1; }
  ENTRY=0x$("$BU/z80-none-elf-nm" bin/DEMO.obj | awk '$3=="__start" {print $1}')
  python3 app_tools/pack_app.py pack bin/body.bin "$ENTRY" "$LA" bin/_DEMO.bin > /dev/null || exit 1
  rm -f bin/AppIns*.8xv bin/AppIns*.bin
  python3 app_tools/make_segments.py bin/_DEMO.bin bin >/dev/null || exit 1
  size=$(stat -c %s bin/_DEMO.bin); segs=$(( (size + 3 + 65231) / 65232 )); last=$(( size + 3 - (segs - 1) * 65232 ))
}
t0=$(date +%s)
link 0
python3 tools/dev/v15_biascheck.py obj/lto.s || { echo REFUSED-bias; exit 1; }
python3 tools/dev/v15_carrycheck.py obj/lto.s || { echo REFUSED-carry; exit 1; }
pack_all
if [ "$last" -lt 1024 ]; then PAD=$((1024 - last + 64)); link 0; pack_all; fi
echo "[build] OK in $(( $(date +%s)-t0 ))s size=$size room=$(( 0x2AFFFD - size ))"
if [ "$size" -gt $((0x2AFFFD)) ]; then echo "[build] TOO BIG: $size B > 43 flash pages: INST will fail"; exit 3; fi
if [ "$last" -lt 1024 ]; then echo "[build] the last AppIns holds $last B (< 1 KB): INST stops on it"; exit 3; fi
