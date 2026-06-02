#!/usr/bin/env bash
# Build libenzomodules_grid (the Session/grid-method bridge) on macOS, linked
# against the locally-built libenzo_p8_b8.dylib. Darwin counterpart of
# build_grid.sh (which targets Ubuntu). Run with bash (NOT zsh — needs word-split).
#   bash EnzoModules/deps/build_grid_darwin.sh
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENZO="$repo/src/enzo"
GFLIB="$(dirname "$(gfortran -print-file-name=libgfortran.dylib)")"
HDF5=/opt/homebrew
OUT="$repo/EnzoModules/deps/libenzomodules_grid.dylib"

# ---- stage 1: full Enzo serial shared library (slow; skipped if present) -----
# macOS specifics vs the stock/Ubuntu build: real Homebrew gcc-15 (Apple clang
# rejects Enzo's "%"ISYM literals), -ffixed-line-length-132 (fixed-form Fortran),
# and Homebrew HDF5 (arm64; its H5version.h still provides the v16 API that
# -DH5_USE_16_API needs).
libenzo="$(ls "$ENZO"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
if [ -z "$libenzo" ]; then
  echo "[grid] building full Enzo serial library (slow)..."
  FF="-fallow-argument-mismatch -fno-second-underscore -m64 -ffixed-line-length-132"
  ( cd "$repo" && ./configure )
  ( cd "$ENZO"
    make machine-darwin
    make use-mpi-no precision-64 particles-64 integers-32
    make clean
    LIBRARY_PATH="$GFLIB:${LIBRARY_PATH:-}" make lib -j"$(sysctl -n hw.ncpu)" \
      MACH_CXX_NOMPI=g++-15 MACH_CC_NOMPI=gcc-15 MACH_LD_NOMPI=g++-15 \
      MACH_FFLAGS="$FF" MACH_F90FLAGS="$FF" \
      LOCAL_HDF5_INSTALL="$HDF5" LOCAL_FC_INSTALL="$GFLIB" \
      MACH_SHARED_FLAGS=-fPIC SHARED_OPT=-shared )
  libenzo="$(ls "$ENZO"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
fi
[ -n "$libenzo" ] || { echo "ERROR: libenzo not built ($ENZO/libenzo_p*_b*.dylib)"; exit 1; }
libname="$(basename "$libenzo" .dylib)"; libname="${libname#lib}"
echo "[grid] using $libenzo  (-l$libname)"

# Enzo's exact -D defines (so the bridge is ABI-compatible with libenzo).
DFLAGS="$(cd "$ENZO" && make -n Grid_EnzoModulesFixture.o MACH_CXX_NOMPI=g++-15 \
          LOCAL_HDF5_INSTALL="$HDF5" 2>/dev/null \
          | grep -oE '\-D[A-Za-z0-9_=]+' | sort -u | tr '\n' ' ')"
INC="-I$ENZO -I$ENZO/hydro_rk -I$HDF5/include"

objs=()
for src in Grid_EnzoModulesFixture enzomodules_grid_bridge enzomodules_problem_bridge \
           enzomodules_chemistry_bridge enzomodules_radiation_bridge enzomodules_ppm_grid_bridge \
           enzomodules_timing_init enzomodules_amr_bridge enzomodules_mhdct_bridge \
           enzomodules_hierarchy_bridge enzomodules_halo_bridge; do
  from="$repo/EnzoModules/src"; [ -f "$ENZO/$src.C" ] && from="$ENZO"
  echo "[grid] CXX $src"
  g++-15 $DFLAGS $INC -fPIC -O2 -c "$from/$src.C" -o "/tmp/em_$src.o"
  objs+=("/tmp/em_$src.o")
done

echo "[grid] LINK $OUT"
g++-15 -dynamiclib -fPIC -o "$OUT" "${objs[@]}" \
  -L"$ENZO" "-l$libname" -L"$HDF5/lib" -lhdf5 -L"$GFLIB" -lgfortran -lstdc++ \
  -Wl,-rpath,"$ENZO" -Wl,-rpath,"$HDF5/lib" -Wl,-rpath,"$GFLIB"
rm -f "${objs[@]}"
# libenzo's install-name is a bare filename; rewrite the dependency to its
# absolute path so the bridge dylib loads without DYLD_LIBRARY_PATH.
install_name_tool -change "$(basename "$libenzo")" "$libenzo" "$OUT"
echo "[grid] OK -> $OUT"
