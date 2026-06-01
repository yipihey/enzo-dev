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
    hydro.py              #   ergonomic kernel wrappers (twoshock, ...)
  deps/build_pilot.sh     # builds libenzomodules_pilot.so (bridge + leaf kernels)
  fixtures/Hydro/twoshock/*.fixture   # golden inputs/outputs (committed)
  tools/capture_twoshock.py           # capture + validate fixtures
  tests/                  # pytest replay suite
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

The pilot build compiles **only** the bridge plus the leaf Fortran kernels it
exposes — not the full Enzo executable — so it is fast and hermetic.

## The per-kernel workflow

Each kernel goes through the same pipeline; `twoshock` (two-shock approximate
Riemann solver) is the worked pilot:

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
