#!/usr/bin/env bash
#
# Build the EnzoModules pilot shared library.
#
# This compiles ONLY the bridge shim plus the handful of leaf Fortran
# kernels it exposes -- it deliberately does NOT build the full Enzo
# executable (no HDF5/MPI/~1000 C++ files).  That keeps the pilot fast,
# hermetic, and easy to reproduce while we prove out the wrap/test/capture
# pipeline.  Full integration into Enzo's `make lib` target is a follow-up
# once more kernels (and the C++ grid-fixture shims) are wrapped.
#
# Precision contract: double baryons (CONFIG_BFLOAT_8) + 32-bit Fortran
# ints (SMALL_INTS).  EnzoModules queries this at load time via
# enzomodules_*_precision_bytes() and refuses to run on a mismatch.
#
# Usage:   deps/build_pilot.sh [output_so_path]
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
enzo_src="${ENZO_SRC:-$(cd "${here}/../../src/enzo" && pwd)}"
out="${1:-${here}/libenzomodules_pilot.so}"
work="$(mktemp -d)"
trap 'rm -rf "${work}"' EXIT

FC="${FC:-gfortran}"
CXX="${CXX:-g++}"

# Match Enzo's precision configuration via cpp defines.
DEFS="-DCONFIG_BFLOAT_8 -DCONFIG_PFLOAT_8 -DSMALL_INTS"
FFLAGS="-cpp ${DEFS} -I${enzo_src} -O2 -fPIC"
CXXFLAGS="-O2 -fPIC -I${enzo_src}"

# Leaf Fortran kernels exposed by the bridge.  Add files here as more
# kernels are wrapped (flux_twoshock.F, intvar.F, ...).
FKERNELS=(
  "twoshock.F"
)

echo "[build_pilot] ENZO_SRC = ${enzo_src}"
echo "[build_pilot] output   = ${out}"

objs=()
for f in "${FKERNELS[@]}"; do
  obj="${work}/$(basename "${f%.F}").o"
  echo "[build_pilot] FC  ${f}"
  "${FC}" ${FFLAGS} -c "${enzo_src}/${f}" -o "${obj}"
  objs+=("${obj}")
done

echo "[build_pilot] CXX enzomodules_bridge.C"
"${CXX}" ${CXXFLAGS} -c "${enzo_src}/enzomodules_bridge.C" -o "${work}/bridge.o"
objs+=("${work}/bridge.o")

echo "[build_pilot] LD  ${out}"
"${CXX}" -shared -fPIC -o "${out}" "${objs[@]}" -lgfortran

echo "[build_pilot] OK"
