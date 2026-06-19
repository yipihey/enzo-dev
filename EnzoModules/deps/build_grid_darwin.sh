#!/usr/bin/env bash
# Build libenzomodules_grid (the Session/grid-method bridge) on macOS, linked
# against a locally-built libenzo_p8_b8.dylib.  Two flavors:
#   bash EnzoModules/deps/build_grid_darwin.sh           # serial (default)
#   bash EnzoModules/deps/build_grid_darwin.sh mpi       # MPI (MPItrampoline)
#
# The serial flavor is the default and is byte-compatible with the historical
# build.  The mpi flavor builds an Enzo compiled with -DUSE_MPI against the
# MPItrampoline ABI (the project's standard MPI provider: its mpicc/mpicxx must
# be on PATH; the backend compiler is gcc-15 via MPITRAMPOLINE_{CC,CXX}).  The
# two flavors coexist: serial libenzo lives in src/enzo/, MPI libenzo in
# src/enzo/mpi/, and the bridges are libenzomodules_grid{,_mpi}.dylib.
# Run with bash (NOT zsh — needs word-split).
set -euo pipefail

FLAVOR="${1:-serial}"
case "$FLAVOR" in
  serial|mpi|f32) ;;
  *) echo "usage: $0 [serial|mpi|f32]"; exit 1 ;;
esac

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ENZO="$repo/src/enzo"
GFLIB="$(dirname "$(gfortran -print-file-name=libgfortran.dylib)")"
HDF5=/opt/homebrew
FF="-fallow-argument-mismatch -fno-second-underscore -m64 -ffixed-line-length-132"
JULIA="${JULIA:-$(command -v julia || echo julia)}"

# ---- flavor-specific knobs ---------------------------------------------------
if [ "$FLAVOR" = "mpi" ]; then
  # MPItrampoline is the standard MPI provider.  We compile with gcc-15 directly
  # against the trampoline ABI (its own mpicc wrapper bakes in build-sandbox
  # paths, so we don't use it): -I$TRAMP/include and -L$TRAMP/lib -lmpitrampoline.
  # TRAMP = the MPItrampoline_jll artifact dir (set MPITRAMPOLINE_DIR, or resolve
  # via the EnzoLib test project's MPIPreferences/MPItrampoline_jll).
  TRAMP="${MPITRAMPOLINE_DIR:-}"
  if [ -z "$TRAMP" ]; then
    TRAMP="$("$JULIA" --project="$repo/../Vespa.jl/lib/EnzoLib/test" \
             -e 'import MPItrampoline_jll as T; T.is_available() && print(T.artifact_dir)' 2>/dev/null || true)"
  fi
  [ -n "$TRAMP" ] && [ -f "$TRAMP/include/mpi.h" ] || \
    { echo "ERROR: MPItrampoline artifact not found (set MPITRAMPOLINE_DIR, or configure MPI.jl: MPIPreferences.use_jll_binary(\"MPItrampoline_jll\"))"; exit 1; }
  echo "[grid] MPItrampoline: $TRAMP"
  # IMPORTANT: use the HDF5 keg include (no mpi.h) instead of the umbrella
  # /opt/homebrew/include, whose mpi.h is open-mpi's and would shadow the
  # MPItrampoline header → undefined `ompi_*` symbols at link.
  HDF5=/opt/homebrew/opt/hdf5
  CXX="g++-15"; CC="gcc-15"
  ENZO_LIBDIR="$ENZO/mpi"
  OUT="$repo/EnzoModules/deps/libenzomodules_grid_mpi.dylib"
  MPI_MAKE_TARGET="use-mpi-yes"
  EXTRA_CONFIG=""    # keep TRANSFER on (parity with serial; fixtures need it)
  # -fpermissive: the FLD/FSProb radiation solvers pass int where MPItrampoline's
  # strict MPI_Comm (a pointer) is required — open-mpi's looser typedef hid it.
  # We don't run those solvers; downgrade the conversion to a warning so the
  # library still builds with full feature parity.  Core MPI code is strict-clean.
  #
  # CRITICAL — embedding: compile against the MPItrampoline headers (-I) but DO NOT
  # link libmpitrampoline.  The MPI_* symbols are left undefined and resolved at
  # load time from the host process's already-loaded MPItrampoline (MPI.jl's).
  # Two reasons NOT to link it / NOT to use a blanket `-undefined dynamic_lookup`:
  #   (1) linking libmpitrampoline (two-level) loads a 2nd trampoline instance in
  #       the MPI.jl process → MPItrampoline aborts (forbids double-load);
  #   (2) a blanket dynamic_lookup also floats libenzo's C++ stdlib symbols, which
  #       then mis-resolve across C++ runtimes (gcc-15 libstdc++ vs the macOS libc++
  #       stub) → std::locale double-free abort in static init at dlopen.
  # So we keep two-level namespace (correct C++) and float ONLY the MPI symbols via
  # explicit `-Wl,-U` (the exact set Enzo references; a new MPI call would fail the
  # link with a clear "undefined _MPI_X", prompting an addition here).
  MPI_USYMS=""
  for s in MPI_Abort MPI_Allgather MPI_Allgatherv MPI_Allreduce MPI_Alltoall \
           MPI_Alltoallv MPI_Barrier MPI_Bcast MPI_BYTE MPI_Cancel MPI_CHAR \
           MPI_Comm_create_errhandler MPI_Comm_rank MPI_Comm_set_errhandler \
           MPI_Comm_size MPI_COMM_WORLD MPI_DOUBLE MPI_ERR_OTHER MPI_Errhandler_free \
           MPI_Error_class MPI_Error_string MPI_Finalize MPI_Finalized MPI_FLOAT \
           MPI_Gather MPI_Init MPI_Initialized MPI_INT MPI_Irecv MPI_Isend \
           MPI_LONG_DOUBLE MPI_LONG_INT MPI_LONG_LONG_INT MPI_MAX MPI_MIN MPI_PACKED \
           MPI_PROC_NULL MPI_Recv MPI_Reduce MPI_REQUEST_NULL MPI_Send MPI_Ssend \
           MPI_STATUS_IGNORE MPI_SUCCESS MPI_SUM MPI_Test MPI_Testsome MPI_Type_commit \
           MPI_Type_contiguous MPI_Type_size MPI_Wait MPI_Waitall MPI_Waitsome MPI_Wtime; do
    MPI_USYMS+=" -Wl,-U,_$s"
  done
  MACH_OVERRIDES=(MACH_CXX_MPI="$CXX -fpermissive" MACH_CC_MPI="$CC" \
                  MACH_LD_MPI="$CXX$MPI_USYMS" \
                  LOCAL_MPI_INSTALL="$TRAMP" LOCAL_LIBS_MPI="")
  MPI_INC="-I$TRAMP/include"
  MPI_LINK="$MPI_USYMS"
  PREC_TARGETS="precision-64 particles-64"
elif [ "$FLAVOR" = "f32" ]; then
  # 32-bit baryon + 32-bit particle Enzo (p4_b4) — the faithful-precision CPU
  # reference for the EnzoNG f32 GPU kernels.  Coexists with the f64 serial bridge
  # (which the bit-tight f64 oracle tests depend on): its libenzo lives in
  # $ENZO/f32 and its bridge is libenzomodules_grid_f32.dylib, selected at load
  # time via ENV["ENZOMODULES_GRID_LIB"].  The C-ABI is precision-independent
  # (the bridge casts enzo_float<->double element-wise), so no Julia-side change.
  CXX="g++-15"; CC="gcc-15"
  ENZO_LIBDIR="$ENZO/f32"
  OUT="$repo/EnzoModules/deps/libenzomodules_grid_f32.dylib"
  MPI_MAKE_TARGET="use-mpi-no"
  EXTRA_CONFIG=""
  MACH_OVERRIDES=(MACH_CXX_NOMPI="$CXX" MACH_CC_NOMPI="$CC" MACH_LD_NOMPI="$CXX")
  MPI_INC=""; MPI_LINK=""
  PREC_TARGETS="precision-32 particles-32"
else
  CXX="g++-15"; CC="gcc-15"
  ENZO_LIBDIR="$ENZO"
  OUT="$repo/EnzoModules/deps/libenzomodules_grid.dylib"
  MPI_MAKE_TARGET="use-mpi-no"
  EXTRA_CONFIG=""
  MACH_OVERRIDES=(MACH_CXX_NOMPI="$CXX" MACH_CC_NOMPI="$CC" MACH_LD_NOMPI="$CXX")
  MPI_INC=""; MPI_LINK=""
  PREC_TARGETS="precision-64 particles-64"
fi
mkdir -p "$ENZO_LIBDIR"

# ---- stage 1: full Enzo shared library (slow; skipped if present) ------------
# macOS specifics: real Homebrew gcc-15 (Apple clang rejects Enzo's "%"ISYM
# literals), -ffixed-line-length-132 (fixed-form Fortran), Homebrew HDF5 (arm64;
# its H5version.h still provides the v16 API that -DH5_USE_16_API needs).
libenzo="$(ls "$ENZO_LIBDIR"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
if [ -z "$libenzo" ]; then
  echo "[grid] building full Enzo $FLAVOR library (slow)..."
  ( cd "$repo" && ./configure )
  # The serial dylib (if any) lives at $ENZO/libenzo_*.dylib; a mpi build would
  # clobber that name during `make lib`, so stash and restore it.
  # `make lib` always emits libenzo_*.dylib into $ENZO; a flavor that keeps its
  # library elsewhere ($ENZO/mpi or $ENZO/f32) must stash any serial libenzo there
  # first (so it isn't clobbered) and move the freshly-built one to its own dir.
  stash=""
  if [ "$ENZO_LIBDIR" != "$ENZO" ]; then
    existing="$(ls "$ENZO"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
    if [ -n "$existing" ]; then stash="/tmp/$(basename "$existing").serial-stash"; mv "$existing" "$stash"; fi
  fi
  ( cd "$ENZO"
    make machine-darwin
    make "$MPI_MAKE_TARGET" $PREC_TARGETS integers-32 ${EXTRA_CONFIG:+$EXTRA_CONFIG}
    make grackle-yes      # use_grackle cooling (the high-z chemistry runs need it)
    make clean
    LIBRARY_PATH="$GFLIB:${LIBRARY_PATH:-}" make lib -j"$(sysctl -n hw.ncpu)" \
      "${MACH_OVERRIDES[@]}" \
      MACH_FFLAGS="$FF" MACH_F90FLAGS="$FF" \
      LOCAL_HDF5_INSTALL="$HDF5" LOCAL_FC_INSTALL="$GFLIB" \
      LOCAL_GRACKLE_INSTALL="${GRACKLE_INSTALL:-$HOME/grackle_install_f32}" \
      MACH_SHARED_FLAGS=-fPIC SHARED_OPT=-shared )
  built="$(ls "$ENZO"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
  if [ "$ENZO_LIBDIR" != "$ENZO" ] && [ -n "$built" ]; then mv "$built" "$ENZO_LIBDIR/"; fi
  [ -n "$stash" ] && mv "$stash" "$ENZO/$(basename "$stash" .serial-stash)"
  libenzo="$(ls "$ENZO_LIBDIR"/libenzo_p*_b*.dylib 2>/dev/null | head -1 || true)"
fi
[ -n "$libenzo" ] || { echo "ERROR: libenzo not built ($ENZO_LIBDIR/libenzo_p*_b*.dylib)"; exit 1; }
libname="$(basename "$libenzo" .dylib)"; libname="${libname#lib}"
echo "[grid] using $libenzo  (-l$libname)"

# Pin the Make.config MPI *and precision* state to THIS flavor before extracting
# DFLAGS.  Stage 1 sets these only when it builds libenzo; if libenzo is cached
# (the common case) the config still reflects whatever flavor was built LAST.
# Building the serial bridge while the config says use-mpi-yes would compile the
# bridge objects with -DUSE_MPI and fail the (MPI-less) serial link on _MPI_*; and
# building the f32 bridge while the config says precision-64 would bake the wrong
# enzo_float size into the bridge -D flags (ABI mismatch vs the p4_b4 libenzo) —
# so set both explicitly here.
( cd "$ENZO" && make "$MPI_MAKE_TARGET" $PREC_TARGETS integers-32 >/dev/null 2>&1 || true )

# Enzo's exact -D defines (so the bridge is ABI-compatible with libenzo).  The
# mpi flavor's `make -n` runs under use-mpi-yes, so it picks up -DUSE_MPI and the
# MPI include path automatically.
DFLAGS="$(cd "$ENZO" && make -n Grid_EnzoModulesFixture.o "${MACH_OVERRIDES[@]}" \
          LOCAL_HDF5_INSTALL="$HDF5" 2>/dev/null \
          | grep -oE '\-D[A-Za-z0-9_=]+' | sort -u | tr '\n' ' ')"
INC="-I$ENZO -I$ENZO/hydro_rk -I$HDF5/include"

objs=()
for src in Grid_EnzoModulesFixture enzomodules_grid_bridge enzomodules_problem_bridge \
           enzomodules_chemistry_bridge enzomodules_radiation_bridge enzomodules_ppm_grid_bridge \
           enzomodules_timing_init enzomodules_amr_bridge enzomodules_mhdct_bridge \
           enzomodules_hierarchy_bridge enzomodules_halo_bridge enzomodules_hydro_rk_bridge \
           enzomodules_mg_bridge; do
  from="$repo/EnzoModules/src"; [ -f "$ENZO/$src.C" ] && from="$ENZO"
  echo "[grid] CXX $src"
  $CXX $DFLAGS $INC $MPI_INC -fPIC -O2 -c "$from/$src.C" -o "/tmp/em_$src.o"
  objs+=("/tmp/em_$src.o")
done

echo "[grid] LINK $OUT"
$CXX -dynamiclib -fPIC -o "$OUT" "${objs[@]}" \
  -L"$ENZO_LIBDIR" "-l$libname" -L"$HDF5/lib" -lhdf5 -L"$GFLIB" -lgfortran -lstdc++ \
  $MPI_LINK \
  -Wl,-headerpad_max_install_names \
  -Wl,-rpath,"$ENZO_LIBDIR" -Wl,-rpath,"$HDF5/lib" -Wl,-rpath,"$GFLIB"
rm -f "${objs[@]}"
# libenzo's install-name is a bare filename; rewrite the dependency to its
# absolute path so the bridge dylib loads without DYLD_LIBRARY_PATH.
install_name_tool -change "$(basename "$libenzo")" "$libenzo" "$OUT"
# Note: the mpi flavor intentionally has NO libmpitrampoline dependency — MPI_*
# symbols are resolved at load time from the host's MPItrampoline (see above).
echo "[grid] OK -> $OUT"

# ---- stage 3: the ADR-0005 #3 subprocess worker -----------------------------
# A standalone (non-Julia) host process that dlopens the bridge above and serves
# the EnzoNG bridge over the rpc.jl wire protocol (control channel + shared file).
# Its per-symbol typed dispatch is GENERATED from the bridge manifest so it stays
# in lockstep with session.jl; the contract hash baked into the generated header
# is checked at the handshake.  Carrying no Julia runtime, its gcc/libstdc++ stack
# never meets Julia's libc++ — the collision that blocks in-process MPI cannot
# occur, so this is the process that will host the MPI libenzo (mpiexec -n N).
WORKER_SRC="$repo/EnzoModules/src/enzomodules_worker.C"
WORKER_INC="$repo/EnzoModules/src/enzomodules_worker_dispatch.inc"
# The worker is OPTIONAL (the bridge dylib above is the primary artifact); its
# dispatch is generated by Julia, so skip gracefully if julia isn't resolvable
# (e.g. juliaup not on PATH) rather than failing the whole build under `set -e`.
if [ -f "$WORKER_SRC" ] && "$JULIA" --version >/dev/null 2>&1; then
  echo "[worker] generating dispatch from the bridge manifest"
  "$JULIA" --project="$repo/../Vespa.jl/lib/EnzoLib/test" \
    "$repo/EnzoModules/tools/gen_worker_dispatch.jl" "$WORKER_INC"
  if [ "$FLAVOR" = "mpi" ]; then
    # MPI worker: owns MPI_Init in its own (Julia-free) process, so it links the
    # MPItrampoline directly (no double-load hazard — that was specific to the
    # MPI.jl-in-process case).  Launched `mpiexec -n N`; rank 0 owns the control
    # channel and broadcasts each command so collective bridge calls stay in lockstep.
    WORKER_OUT="$repo/EnzoModules/deps/enzomodules_worker_mpi"
    echo "[worker] CXX enzomodules_worker_mpi (-DUSE_MPI, MPItrampoline)"
    "$CXX" -std=c++17 -O2 -Wall -DUSE_MPI -I"$TRAMP/include" -I"$repo/EnzoModules/src" \
      "$WORKER_SRC" -o "$WORKER_OUT" -L"$TRAMP/lib" -lmpitrampoline -Wl,-rpath,"$TRAMP/lib" -ldl
  else
    WORKER_OUT="$repo/EnzoModules/deps/enzomodules_worker"
    [ "$FLAVOR" = "f32" ] && WORKER_OUT="${WORKER_OUT}_f32"
    echo "[worker] CXX $(basename "$WORKER_OUT") (serial)"
    "$CXX" -std=c++17 -O2 -Wall -I"$repo/EnzoModules/src" \
      "$WORKER_SRC" -o "$WORKER_OUT" -ldl
  fi
  echo "[worker] OK -> $WORKER_OUT"
fi
