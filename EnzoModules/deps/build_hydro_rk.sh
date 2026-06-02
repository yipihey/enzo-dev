#!/usr/bin/env bash
#
# Build the EnzoModules hydro_rk shared library.
#
# Unlike the Fortran-kernel pilot (build_pilot.sh, standalone), the hydro_rk
# C++ solvers depend on Enzo's headers and global state, so this links against
# the full Enzo shared library (libenzo_p8_b8.so).  Two stages:
#
#   1. Build Enzo as a serial shared library (slow; skipped if already built).
#      NOTE: the stock `lib` target only sets shared flags for macOS, so we
#      pass MACH_SHARED_FLAGS=-fPIC and SHARED_OPT=-shared explicitly.
#   2. Compile enzomodules_hydro_rk_bridge.C with Enzo's exact DEFINES and link
#      it against libenzo.
#
# Requires: gfortran, g++, and a serial HDF5 (libhdf5-dev; the ubuntu machine
# config expects /usr/include/hdf5/serial and -lhdf5_serial).
#
# Usage:   deps/build_hydro_rk.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
enzo_src="${repo}/src/enzo"
out="${here}/libenzomodules_hydrork.so"

cd "${repo}"

# ---- stage 1: Enzo serial shared library --------------------------------
libenzo="$(ls "${enzo_src}"/libenzo_p*_b*.so 2>/dev/null | head -1 || true)"
if [ -z "${libenzo}" ]; then
  echo "[build_hydro_rk] building Enzo shared library (this is slow)..."
  ./configure
  ( cd "${enzo_src}"
    make machine-ubuntu
    make use-mpi-no
    make precision-64
    make particles-64
    make integers-32
    make clean
    make lib -j"$(nproc)" MACH_SHARED_FLAGS=-fPIC SHARED_OPT=-shared )
  libenzo="$(ls "${enzo_src}"/libenzo_p*_b*.so 2>/dev/null | head -1)"
fi
[ -n "${libenzo}" ] || { echo "ERROR: Enzo shared library not built"; exit 1; }
libname="$(basename "${libenzo}" .so)"; libname="${libname#lib}"
echo "[build_hydro_rk] using ${libenzo}"

# ---- stage 2: bridge object + link --------------------------------------
# Pull Enzo's exact preprocessor DEFINES so the bridge is ABI-compatible.
defines="$(cd "${enzo_src}" && make -s show-flags 2>/dev/null | sed -n 's/^DEFINES = //p')"
[ -n "${defines}" ] || { echo "ERROR: could not read Enzo DEFINES"; exit 1; }

CXX="${CXX:-g++}"
src_dir="$(cd "${here}/../src" && pwd)"   # our bridge sources (Enzo headers via -I)
echo "[build_hydro_rk] CXX enzomodules_hydro_rk_bridge.C"
${CXX} ${defines} -I"${enzo_src}" -I"${enzo_src}/hydro_rk" \
    -I/usr/include/hdf5/serial -fPIC -O2 \
    -c "${src_dir}/enzomodules_hydro_rk_bridge.C" -o "${here}/hrk_bridge.o"

echo "[build_hydro_rk] LD ${out}"
${CXX} -shared -fPIC -o "${out}" "${here}/hrk_bridge.o" \
    -L"${enzo_src}" -l"${libname}" -lhdf5_serial -lz -lgfortran \
    -Wl,-rpath,"${enzo_src}"
rm -f "${here}/hrk_bridge.o"

echo "[build_hydro_rk] OK -> ${out}"
echo "[build_hydro_rk] (libenzo rpath baked to ${enzo_src})"
