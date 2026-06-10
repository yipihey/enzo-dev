#!/usr/bin/env bash
# Build a standalone shared library exposing Enzo's *actual* solve_rate_cool.F as
# a one-zone, C-callable `enzo_chem_onezone` — the ground-truth reference for the
# ChemistryKernels parity layer. No full Enzo / EnzoLib bridge build required:
# just the chemistry Fortran from src/enzo + the wrapper here, via gfortran.
#
# Usage:  FC=gfortran ./build_enzo_chem.sh [src/enzo dir] [output dir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="${1:-$HERE/../../../../src/enzo}"
OUT="${2:-$HERE/build}"
FC="${FC:-gfortran}"
DEFS="-DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8 -DLARGE_INTS"
FLAGS="-cpp -ffixed-line-length-132 -w -fPIC $DEFS"

mkdir -p "$OUT"; cd "$OUT"
SRCS="calc_rates cool1d_multi solve_rate_cool coll_rates cie_thin_cooling_rate colh2diss interpolate"
for f in $SRCS; do
  $FC $FLAGS -I"$SRC" -c "$SRC/$f.F" -o "$f.o"
done
$FC $FLAGS -I"$SRC" -c "$HERE/enzo_chem_wrap.F" -o enzo_chem_wrap.o
$FC -shared -fPIC -o libenzochem.so enzo_chem_wrap.o \
    calc_rates.o cool1d_multi.o solve_rate_cool.o coll_rates.o \
    cie_thin_cooling_rate.o colh2diss.o interpolate.o
echo "built $OUT/libenzochem.so"
