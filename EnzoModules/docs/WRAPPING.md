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

### 1. Add the C-ABI shim — `src/enzo/enzomodules_bridge.{h,C}`

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

## Notes for composite / C++-driver kernels

The leaf Fortran kernels are pure (array in/out, no globals). The C++ layers
(`Grid_*EulerSweep`, `Grid_SolvePPM_DE`) instead operate on a fully built
`grid` object. To wrap those directly you will need a small **fixture-builder
shim** that constructs the minimal grid/global state the routine reads — the
same idea as the libyt `Grid_ConvertToLibyt` conversion. `ppm_sweep_1d` sidesteps
this by reproducing the *numerical core* of the sweep at the Fortran level; that
is usually the cheaper and more faithful first step.
