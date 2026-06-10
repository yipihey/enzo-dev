#!/usr/bin/env bash
# Build the upstream Grackle submodule into a shared library for the Grackle.jl
# binding/oracle. Requires CMake, a C/Fortran toolchain and HDF5.
#
# Usage:  CC=gcc FC=gfortran HDF5_ROOT=/path ./build_grackle.sh
# Result: deps/build/libgrackle.so  (point Grackle.jl at it, or set GRACKLE_LIB)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/../extern/grackle"
OUT="$HERE/build"
[ -f "$SRC/CMakeLists.txt" ] || { echo "init the submodule: git submodule update --init $SRC"; exit 1; }

cmake -B "$OUT/_cmake" -S "$SRC" \
  -DCMAKE_INSTALL_PREFIX="$OUT/_install" \
  -DGRACKLE_USE_DOUBLE=ON -DBUILD_SHARED_LIBS=ON \
  ${HDF5_ROOT:+-DHDF5_ROOT="$HDF5_ROOT"} \
  ${CC:+-DCMAKE_C_COMPILER="$CC"} ${FC:+-DCMAKE_Fortran_COMPILER="$FC"}
cmake --build "$OUT/_cmake" -j"${JOBS:-4}"
cmake --install "$OUT/_cmake"
# expose a stable libgrackle.so symlink for the binding
lib=$(find "$OUT/_install" -name "libgrackle*.so" | head -1)
ln -sf "$lib" "$OUT/libgrackle.so"
echo "built $OUT/libgrackle.so -> $lib"
