# EnzoModules

Extracted, fixture-certified, **callable legacy Enzo compute kernels**.

EnzoModules builds selected legacy C/Fortran Enzo kernels as a small shared
library, exposes each one as a documented Python function (via `ctypes`),
and ships **golden fixtures** — captured legacy inputs/outputs — plus a
tolerance/diff utility used to certify them.

It exists to make rewriting Enzo (e.g. EnzoNG, in a separate repository)
**safe and incremental**: you can call the reference legacy kernel directly,
and verify any rewrite against the captured outputs with a shared tolerance
policy. The two durable artifacts — the compiled bridge `.so` and the
plain-text fixtures — are **language-neutral**, so a Python *or* a Julia
rewrite can reuse them.

> Scope: this package only *extracts, wraps, and certifies the legacy
> kernels*. The native rewrites, backend dispatch ("keep both, use the
> fastest"), and rewrite-vs-legacy benchmarking live in the rewrite repo and
> reuse the fixtures + `diff` from here.

## Layout

```
EnzoModules/
  enzomodules/            # the importable package
    bridge.py             #   ctypes loader + precision check + raw bindings
    diff.py               #   Tolerance / isclose / compare
    fixtures.py           #   read/write the .fixture format
    hydro.py              #   ergonomic kernel wrappers (twoshock, ppm_sweep_1d, ...)
    examples/             #   worked examples built on the certified kernels
      riemann.py          #     exact Riemann solver (analytic truth)
      ppm_sod.py          #     full 1D PPM hydro driver (Sod shock tube)
  deps/build_pilot.sh     # builds libenzomodules_pilot.so (bridge + kernel closure)
  fixtures/Hydro/*/*.fixture          # golden inputs/outputs (committed)
  tools/capture_*.py      # capture + validate fixtures
  tools/demo_sod.py       # run the Sod tube and compare to exact
  tests/                  # pytest replay + physics suite
  docs/WRAPPING.md        # how to wrap & certify a new kernel (start here to contribute)
```

The C side is two files in the Enzo tree:
`src/enzo/enzomodules_bridge.{h,C}` — thin `extern "C"` shims that forward to
the existing Fortran kernels (nothing numerical is reimplemented). This
follows the libyt integration pattern (`ExposeHierarchyToLibyt.C`): a narrow
C ABI over the in-tree implementation.

## Quick start

```bash
cd EnzoModules
./deps/build_pilot.sh            # needs gfortran + g++ only (no HDF5/MPI)
python tools/capture_twoshock.py # run kernel, validate vs exact solver, write fixtures
python -m pytest -q              # replay tests
```

The pilot build compiles **only** the bridge plus the Fortran kernels it
exposes (and their internal call-closure) — not the full Enzo executable — so
it is fast and hermetic.

See the full PPM solver in action against the analytic solution:

```bash
python tools/demo_sod.py            # ASCII plot: PPM vs exact Riemann, L1 error
```

## Two worked examples

EnzoModules ships two reference wrappings (see `docs/WRAPPING.md` for the
step-by-step recipe to add your own):

- **`twoshock`** — a *leaf* Fortran kernel (the two-shock approximate Riemann
  solver): pure array-in/array-out, the simplest case.
- **`ppm_sweep_1d`** — a *composite* solver reproducing the numerical core of
  Enzo's Eulerian PPM directional sweep (`Grid_xEulerSweep.C`):
  `inteuler → twoshock → flux_twoshock → euler`. The
  `enzomodules.examples.ppm_sod` driver evolves a Sod shock tube using only
  this one certified call and matches the exact Riemann solution to an L1
  density error of ~1.5e-3 on 200 cells.

## C++ solvers: hydro_rk and Dedner MHD

The Fortran-kernel examples above build a tiny standalone library.  The Enzo
C++ solvers (`hydro_rk`, ZEUS) instead depend on Enzo's headers and global
state, so they are wrapped by linking against the **full Enzo shared
library**.  A second bridge (`src/enzo/enzomodules_hydro_rk_bridge.C`) and
library (`libenzomodules_hydrork.so`) cover this:

```bash
./deps/build_hydro_rk.sh    # builds Enzo as a serial .so (slow), then the bridge
./deps/build_grid.sh        # grid-method bridge (reuses the Enzo .so)
python -m pytest tests/test_hydro_rk.py tests/test_zeus.py -q
```

Wrapped and certified so far:

- **hydro_rk Riemann solvers** (`enzomodules.bridge.hydro_rk_line`): the
  HLL / HLLC / LLF + PLM 1D line solvers.  `examples.hydro_rk_sod` evolves a
  Sod tube with each and matches the exact Riemann solution to L1(rho) ~
  2.5–3.2e-3 on 200 cells.
- **Dedner divergence-cleaning MHD** (`enzomodules.bridge.mhd_rk_line`): the
  HLLD (and HLL/LLF) + PLM MHD line solver with GLM cleaning.
  `examples.mhd_brio_wu` evolves the Brio & Wu shock tube and reproduces its
  structure while keeping the normal field `Bx` constant (divergence-free).
- **ZEUS** (`enzomodules.bridge.zeus_sweep_1d`): the operator-split
  finite-difference `grid::ZeusSolver`.  `examples.zeus_sod` evolves a Sod tube
  and matches the exact solution (a bit more diffusive, as expected from ZEUS
  artificial viscosity).

Two build fixes make this possible (both in `deps/`): `MACH_SHARED_FLAGS`/
`SHARED_OPT` for a non-macOS Enzo shared build, and the grid-fixture methods on
the `grid` class (below).

### The grid-fixture infrastructure

ZEUS is a `grid::` *method* — it reads a fully constructed `grid` object, not
plain arrays.  Wrapping it (and any other grid-method solver) is done with a
small, generic set of methods added to the `grid` class
(`src/enzo/Grid_EnzoModulesFixture.C`, declared in `Grid.h` next to the libyt
hooks):

| method | purpose |
|---|---|
| `EnzoModulesSetupGrid(rank, dims, left, right, nfields, types, dt)` | dimensions, fields, allocation, timestep, mark local |
| `EnzoModulesSetField` / `EnzoModulesGetField` | copy a field in / out |
| `EnzoModulesFieldIndex(type)` | locate a field by `FieldType` |

The grid bridge (`src/enzo/enzomodules_grid_bridge.C`) seeds globals with
Enzo's own `SetDefaultGlobalValues`, builds a fixture, calls the method, and
reads the result back.  **This is the foundation for wrapping radiation
transfer, gravity/Poisson, chemistry, and the rest of the grid-method
solvers** — they reuse the same primitives; see `docs/WRAPPING.md`.

> Still to do: the multi-dimensional sweeps (y/z), and applying the grid
> fixture to radiation transfer and gravity (the next solver families).

## The per-kernel workflow

**To add a kernel, follow [`docs/WRAPPING.md`](docs/WRAPPING.md)** — the
step-by-step recipe, with `twoshock` (leaf), `ppm_sweep_1d` (composite Fortran)
and `hydro_rk_line` (libenzo-linked C++) as worked references.

In short, each kernel goes through the same pipeline:

1. **Wrap** — add an `extern "C"` shim in `enzomodules_bridge.{h,C}` that
   forwards to the Fortran routine; add the `.F` file to `build_pilot.sh`.
2. **Bind** — add a `ctypes` signature in `bridge.py` and an ergonomic
   wrapper in the relevant submodule (`hydro.py`, ...).
3. **Capture** — drive the kernel on canonical inputs, validate the result
   against an independent reference (here, an exact Riemann solver), and
   serialize `(inputs, params, outputs)` to `fixtures/.../*.fixture`.
4. **Test** — replay fixtures through the wrapped kernel and assert it
   reproduces the reference within the kernel's tolerance; include a negative
   control.

Then climb the call graph: `flux_twoshock`, `intvar`, `inteuler`, `intprim`,
`euler` → the C++ `Grid_{x,y,z}EulerSweep` → `Grid_SolvePPM_DE`. The upper
(C++) layers will need a small grid/global *fixture-construction* shim,
modelled again on the libyt `ConvertToLibyt` conversion.

## Precision contract

The pilot library is built with **double** baryon fields (`CONFIG_BFLOAT_8`)
and **32-bit** Fortran integers (`SMALL_INTS`). The library reports its own
precision (`enzomodules_*_precision_bytes`), and `bridge.check_precision()`
refuses to run on a mismatch — precision drift is exactly the class of silent
bug this framework is meant to surface, not hide.

## Reusing fixtures from another language

The `.fixture` format is line-oriented text (`key = value` scalars and
`@name v…` float arrays). A Julia (or any) rewrite can read these directly to
verify its own port without depending on Python — it only needs the same
tolerance policy (mirror `diff.py`).

## CI

`.github/workflows/enzomodules.yml` builds the pilot library, regenerates and
validates the fixtures (failing on drift from the analytic solution or from
the committed fixtures), and runs the replay tests. It is independent of the
CircleCI regression suite so a slow full-Enzo build never gates this fast
feedback.
