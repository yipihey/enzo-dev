#!/usr/bin/env bash
#
# Build the EnzoModules grid-method shared library (libenzomodules_grid.so).
#
# Wraps grid:: solver *methods* (ZEUS, and later radiation transfer, gravity)
# via the generic grid fixture primitives.  Like build_hydro_rk.sh it links
# against the full Enzo shared library; it additionally compiles the grid
# fixture methods (Grid_EnzoModulesFixture.C).  Adding those non-virtual
# methods to Grid.h does not change object layout, so linking against an
# already-built libenzo is ABI-safe.
#
# Usage:   deps/build_grid.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../.." && pwd)"
enzo_src="${repo}/src/enzo"
out="${here}/libenzomodules_grid.so"

cd "${repo}"

# ---- Enzo serial shared library (shared with build_hydro_rk.sh) ----------
libenzo="$(ls "${enzo_src}"/libenzo_p*_b*.so 2>/dev/null | head -1 || true)"
if [ -z "${libenzo}" ]; then
  echo "[build_grid] building Enzo shared library (slow)..."
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
echo "[build_grid] using ${libenzo}"

defines="$(cd "${enzo_src}" && make -s show-flags 2>/dev/null | sed -n 's/^DEFINES = //p')"
[ -n "${defines}" ] || { echo "ERROR: could not read Enzo DEFINES"; exit 1; }

CXX="${CXX:-g++}"
inc="-I${enzo_src} -I${enzo_src}/hydro_rk -I/usr/include/hdf5/serial"

for src in Grid_EnzoModulesFixture enzomodules_grid_bridge enzomodules_problem_bridge enzomodules_chemistry_bridge enzomodules_radiation_bridge enzomodules_ppm_grid_bridge enzomodules_timing_init; do
  echo "[build_grid] CXX ${src}.C"
  ${CXX} ${defines} ${inc} -fPIC -O2 -c "${enzo_src}/${src}.C" -o "${here}/${src}.o"
done

echo "[build_grid] LD ${out}"
${CXX} -shared -fPIC -o "${out}" \
    "${here}/Grid_EnzoModulesFixture.o" "${here}/enzomodules_grid_bridge.o" "${here}/enzomodules_problem_bridge.o" "${here}/enzomodules_chemistry_bridge.o" "${here}/enzomodules_radiation_bridge.o" "${here}/enzomodules_ppm_grid_bridge.o" "${here}/enzomodules_timing_init.o" \
    -L"${enzo_src}" -l"${libname}" -lhdf5_serial -lz -lgfortran \
    -Wl,-rpath,"${enzo_src}"
rm -f "${here}/Grid_EnzoModulesFixture.o" "${here}/enzomodules_grid_bridge.o" "${here}/enzomodules_problem_bridge.o" "${here}/enzomodules_chemistry_bridge.o" "${here}/enzomodules_radiation_bridge.o" "${here}/enzomodules_ppm_grid_bridge.o" "${here}/enzomodules_timing_init.o"

echo "[build_grid] OK -> ${out}"
