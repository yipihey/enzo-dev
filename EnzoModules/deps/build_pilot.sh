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
# Fortran flags mirror Enzo's linux-gnu machine config (Make.mach.linux-gnu):
# fixed-line-length-132 (some kernels emit long preprocessed lines, e.g.
# CALL f_warning(__FILE__,...)), legacy std, and single trailing underscore
# so the FORTRAN_NAME(NAME)=NAME_ convention in the bridge resolves.
FFLAGS="-cpp ${DEFS} -I${enzo_src} -O2 -fPIC -ffixed-line-length-132 -std=legacy -fno-second-underscore"
CXXFLAGS="-O2 -fPIC -I${enzo_src} -DENZOMODULES_STANDALONE"

# Fortran kernels exposed by the bridge, plus their internal call-closure
# (e.g. inteuler -> intvar/intprim/calc_eigen/intpos; flux_twoshock ->
# flux_hll).  Add files here as more kernels are wrapped.
FKERNELS=(
  "twoshock.F"
  # PPM 1D sweep closure (inteuler -> twoshock -> flux_twoshock -> euler):
  "inteuler.F"
  "intvar.F"
  "intprim.F"
  "calc_eigen.F"
  "intpos.F"
  "flux_twoshock.F"
  "flux_hll.F"
  "euler.F"
)
# Note: the ERROR_MESSAGE/WARNING_MESSAGE macros in the reconstruction
# kernels resolve to fc_error/fc_warning, which the bridge provides as light
# standalone stubs (guarded by -DENZOMODULES_STANDALONE) so we avoid pulling
# in Enzo's MPI-dependent c_message.C.

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
src_dir="$(cd "${here}/../src" && pwd)"   # our bridge source (Enzo headers via -I)
"${CXX}" ${CXXFLAGS} -c "${src_dir}/enzomodules_bridge.C" -o "${work}/bridge.o"
objs+=("${work}/bridge.o")

echo "[build_pilot] LD  ${out}"
"${CXX}" -shared -fPIC -o "${out}" "${objs[@]}" -lgfortran

echo "[build_pilot] OK"
