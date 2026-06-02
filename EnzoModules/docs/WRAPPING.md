# How to wrap and certify a legacy Enzo kernel

This is the repeatable recipe for bringing a legacy C/Fortran Enzo routine
into EnzoModules so it can be called, certified against golden fixtures, and
reused by a rewrite. Two worked examples ship in the tree:

| example | kind | files |
|---|---|---|
| `twoshock` | **leaf** Fortran kernel (pure array in/out) | `twoshock.F` |
| `ppm_sweep_1d` | **composite** solver (a pipeline of kernels) | `inteuler` → `twoshock` → `flux_twoshock` → `euler` |

Follow whichever matches your target. The steps are the same; a composite
just orchestrates several `CALL`s inside one shim, mirroring the C++ driver it
replaces.

## The pattern

```
legacy .F / .C          enzomodules_bridge.{h,C}        enzomodules/*.py
─────────────────       ────────────────────────        ──────────────────
subroutine foo(...)  →  extern "C" enzomodules_foo  →   bridge.foo_raw (ctypes)
                        (forwards via FORTRAN_NAME)      hydro.foo  (ergonomic)
                                                         │
                                          tools/capture_foo.py  →  fixtures/.../*.fixture
                                                         │
                                          tests/test_foo.py  (replay + physics)
```

### 1. Add the C-ABI shim — `EnzoModules/src/enzomodules_bridge.{h,C}`

Declare a Fortran prototype (all args by reference) and a thin `extern "C"`
shim that forwards to it via `EM_FORTRAN_NAME` (the trailing-underscore
mangling that matches `macros_and_parameters.h`). Pass scalars by value in the
public ABI and take their address when calling Fortran. **Reimplement
nothing** — the shim only marshals arguments.

For a *composite* kernel, allocate the scratch slabs the C++ driver allocates
(see `Grid_xEulerSweep.C`) and `CALL` the kernels in order. Keep disabled
features (gravity, dual energy, colours…) allocated-but-inert with their flags
set to zero, exactly as the driver does.

Expose the build precision via the existing `enzomodules_*_precision_bytes`
so callers can detect a mismatch.

### 2. Add the Fortran files to the build — `deps/build_pilot.sh`

Append your `.F` file **and its internal call-closure** to `FKERNELS` (e.g.
`inteuler` pulls in `intvar`, `intprim`, `calc_eigen`, `intpos`;
`flux_twoshock` pulls in `flux_hll`). Find the closure with:

```bash
grep -iE '^\s+call ' yourkernel.F
```

If a kernel calls Enzo's error hooks (`ERROR_MESSAGE` / `WARNING_MESSAGE` →
`f_error`/`f_warning`), they are already provided as standalone stubs guarded
by `-DENZOMODULES_STANDALONE`; you do not need `c_message.C`/MPI.

Build and check the symbol resolves:

```bash
./deps/build_pilot.sh
nm -D deps/libenzomodules_pilot.so | grep enzomodules_yourkernel
```

### 3. Add the Python binding — `enzomodules/bridge.py`

Declare `argtypes`/`restype` for `enzomodules_yourkernel` in `_load()`, then a
`yourkernel_raw(...)` that builds `ctypes` arrays, calls, and returns plain
lists. Add an ergonomic wrapper in the relevant submodule
(`enzomodules/hydro.py`, …) with keyword arguments and sensible defaults.

### 4. Capture + validate fixtures — `tools/capture_yourkernel.py`

Drive the kernel on canonical inputs, **validate against an independent
reference** (an analytic solution, an exact solver — see
`enzomodules/examples/riemann.py`), and serialize inputs+outputs with
`fixtures.save_fixture`. The capture script should *fail* if the kernel drifts
from the reference, so regenerating fixtures is itself a check.

### 5. Add tests — `tests/test_yourkernel.py`

- **Replay**: run the wrapped kernel on each fixture's stored inputs and assert
  it reproduces the stored outputs within a tight `Tolerance` (same library →
  effectively bitwise).
- **Negative control**: corrupt a reference and assert the comparison fails.
- **Physics** (recommended): for composite kernels, an end-to-end check against
  analytic truth (e.g. Sod L1 error), like `test_sod_matches_exact_riemann`.

Guard the module with
`pytest.mark.skipif(not bridge.available(), ...)` so a checkout without the C
build still passes.

## Choosing tolerances

`enzomodules/diff.py` mirrors Enzo's yt answer-testing philosophy
(relative + absolute). Use **bitwise** (`BITWISE`) for integer/index logic and
for replay against the same library; use a small `rtol` when comparing across
compilers/precisions or against an analytic solution. State the tolerance and
why in the test.

## Wrapping C++ solvers (link against the full Enzo library)

Fortran leaf kernels are pure and build into a tiny standalone library. The
Enzo C++ solvers (`hydro_rk`, ZEUS) `#include "Grid.h"` + `global_data.h` and
read Enzo's global state, so they are wrapped differently: compile a bridge TU
with Enzo's headers and **link it against the full Enzo shared library**, which
supplies every global's definition. `enzomodules_hydro_rk_bridge.C` +
`deps/build_hydro_rk.sh` are the worked example.

Recipe for a **free-function** C++ solver (e.g. the `hydro_rk` Riemann solvers):

1. Build Enzo as a serial shared library. The stock `lib` target only sets
   shared flags on macOS, so pass `MACH_SHARED_FLAGS=-fPIC SHARED_OPT=-shared`
   (see `build_hydro_rk.sh`). This produces `libenzo_p8_b8.so`.
2. In a bridge TU compiled with Enzo's exact `DEFINES` (read them with
   `make -s show-flags`), `#include` the Enzo headers, declare the solver
   prototype, and write an `extern "C"` shim that **sets the globals the solver
   reads** then calls it. For the `hydro_rk` line solvers those globals are
   `RiemannSolver`, `ReconstructionMethod`, `Gamma`, `EOSType`, `Theta_Limiter`,
   `NumberOfGhostZones`, `NEQ_HYDRO`/`NEQ_MHD`, the prim/flux index globals
   (`iden`, `ivx`, …, `iD`, `iS1`, …), `SmallRho`/`SmallP`, and — for Dedner
   MHD — the cleaning wave speed `C_h`.
3. Link the bridge object against `libenzo` (`-lenzo_p8_b8 -lhdf5_serial -lz
   -lgfortran`, with an rpath to `src/enzo`).
4. Bind with `ctypes` (`bridge.hydro_rk_line` / `mhd_rk_line`), build a
   finite-volume driver (`examples/hydro_rk_sod.py`, `mhd_brio_wu.py`), and
   validate against analytic truth (exact Riemann; Brio-Wu structure +
   divergence-free `Bx`).

Mind the **data conventions** — they are per-solver. The `hydro_rk` Riemann
solvers take primitives in the order `[rho, eint, vx, vy, vz]` (internal energy,
*not* total), `[+Bx, By, Bz, Phi]` for MHD; the conserved energy is
`rho*(eint + 0.5 v^2) + 0.5 B^2`. Always read the solver body to confirm.

## Wrapping C++ `grid::` methods (ZEUS, radiation, gravity) — the grid fixture

ZEUS (`grid::ZeusSolver`), radiation transfer, gravity, etc. are **member
functions** that read a constructed `grid` (`BaryonField`, `GridDimension`,
`GridStartIndex`, boundary state) plus globals. The generic **grid-fixture
infrastructure** handles the grid construction so you don't have to:

`EnzoModules/src/Grid_EnzoModulesFixture.C` adds these methods to the `grid` class
(declared in `Grid.h` next to the libyt hooks; adding non-virtual methods does
not change object layout, so they link against a prebuilt `libenzo`):

- `EnzoModulesSetupGrid(rank, dims, left, right, nfields, field_types, dt)` —
  set the field list, call `PrepareGrid` + `AllocateGrids`, set `dtFixed`, and
  mark the grid local. Returns the cell count.
- `EnzoModulesSetField(i, data)` / `EnzoModulesGetField(i, data)` — copy a field
  in / out (flat `double` arrays).
- `EnzoModulesFieldIndex(field_type)` — locate a field by `FieldType`.

To wrap a new grid-method solver (e.g. an FLD radiation step or the Poisson
solver), add an `extern "C"` shim to `enzomodules_grid_bridge.C` that:

1. calls `SetDefaultGlobalValues(MetaData)` to seed every global, then overrides
   the few the solver needs (e.g. `RadiativeTransferFLD`, `HydroMethod`, …);
2. `EnzoModulesSetupGrid` with the field list the solver expects (add the right
   `FieldType`s — `RadiationFreq0`, `Potential`, species densities, …);
3. fills inputs via `EnzoModulesSetField`, calls the `grid::` method, reads
   outputs via `EnzoModulesGetField`;
4. is bound in `bridge.py`, driven by an `examples/` script, and certified
   against an analytic/reference solution like the others.

`enzomodules_zeus_sweep_1d` (+ `examples/zeus_sod.py`, `tests/test_zeus.py`) is
the worked example. Multi-dimensional sweeps and the other solver families
follow the same recipe — the fixture, bridge, ctypes, and driver scaffolding
are all in place.
