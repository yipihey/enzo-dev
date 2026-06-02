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
  src/                    # our C-ABI bridge + grid-fixture sources (all "ours")
    enzomodules_bridge.{h,C}        #   extern "C" shims over the Fortran kernels
    enzomodules_*_bridge.C          #   grid / problem / radiation / halo / ... bridges
    Grid_EnzoModulesFixture.C       #   grid-class fixture methods
  enzomodules/            # the importable package
    bridge.py             #   ctypes loader + precision check + raw bindings
    problems.py           #   Session / Problem drivers (full step loop in Python)
    diff.py               #   Tolerance / isclose / compare
    fixtures.py           #   read/write the .fixture format
    hydro.py              #   ergonomic kernel wrappers (twoshock, ppm_sweep_1d, ...)
    examples/             #   worked examples built on the certified kernels
      riemann.py          #     exact Riemann solver (analytic truth)
      ppm_sod.py          #     full 1D PPM hydro driver (Sod shock tube)
  deps/build_*.sh         # build the .so libraries (compile src/ against the Enzo headers)
  fixtures/Hydro/*/*.fixture          # golden inputs/outputs (committed)
  tools/capture_*.py      # capture + validate fixtures
  tools/demo_sod.py       # run the Sod tube and compare to exact
  tests/                  # pytest replay + physics suite
  docs/WRAPPING.md        # how to wrap & certify a new kernel (start here to contribute)
```

The C side lives entirely under `EnzoModules/src/` (kept out of the Enzo source
tree so a normal `make enzo` is untouched and upstream stays pristine).
`EnzoModules/src/enzomodules_bridge.{h,C}` are thin `extern "C"` shims that
forward to the existing Fortran kernels (nothing numerical is reimplemented),
following the libyt integration pattern (`ExposeHierarchyToLibyt.C`): a narrow C
ABI over the in-tree implementation, compiled against Enzo's headers via `-I`.

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
library**.  A second bridge (`EnzoModules/src/enzomodules_hydro_rk_bridge.C`) and
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
(`EnzoModules/src/Grid_EnzoModulesFixture.C`, declared in `Grid.h` next to the libyt
hooks):

| method | purpose |
|---|---|
| `EnzoModulesSetupGrid(rank, dims, left, right, nfields, types, dt)` | dimensions, fields (any `FieldType`, incl. species & radiation), allocation, timestep, mark local |
| `EnzoModulesSetField` / `EnzoModulesGetField` / `EnzoModulesFieldIndex` | copy a field in / out, locate by `FieldType` |
| `EnzoModulesSetupParticles(n, nattr)` + position/velocity/mass setters | a **full grid** with particles |
| `EnzoModulesDepositParticles` / `GetDepositField` | run CIC particle-mesh deposit and read it back |

The grid bridge (`EnzoModules/src/enzomodules_grid_bridge.C`) seeds globals with
Enzo's own `SetDefaultGlobalValues`, builds a fixture, calls the method, and
reads the result back.  Certified on this infrastructure so far:

- **ZEUS** (`bridge.zeus_sweep_1d`) — `grid::ZeusSolver`, Sod vs exact.
- **Particles / CIC deposit** (`bridge.cic_deposit`) —
  `grid::DepositParticlePositions`.  `tests/test_particles.py` certifies the
  cloud-in-cell invariants (a particle equidistant from 8 cells splits 1/8
  each, plus linearity, additivity, and mass conservation independent of
  position).
- **Radiation transport / ray tracing** (`bridge.raytrace_uniform`,
  `bridge.hi_cross_section`) — fires a photon package through a uniform-HI grid
  with the legacy ray-tracer (`grid::WalkPhotonPackage`) and reproduces
  **Beer-Lambert** attenuation `N(L) = N0·exp(-n_HI·σ·L)` to <0.5% across
  optical depths (σ from Enzo's own `FindCrossSection`), with full absorption
  in the optically-thick limit.  Also `bridge.rt_identify` — a grid carrying
  the `kphHI`/`PhotoGamma` rate fields, identified by
  `grid::IdentifyRadiativeTransferFields` (`tests/test_radiation.py`).
- **Gravity / Poisson** (`bridge.poisson_solve`) — Enzo's multigrid solver
  (the engine behind `grid::SolveForPotential`).  `tests/test_gravity.py`
  certifies it against manufactured solutions in 1D/2D/3D (recovers
  `sin(pi x)` shapes to <5e-3, converged residual ~1e-9) plus linearity and
  superposition.
- **Chemistry / cooling** (`bridge.chemistry_step`) — the non-Grackle
  primordial network + radiative cooling (`grid::SolveRateAndCoolEquations`;
  rates computed analytically by `InitializeRateData`, no data files) at **all
  three MultiSpecies levels**: 6 species (H, He), 9 species (+ H-/H2/H2+), and
  12 species (+ D/D+/HD).  `tests/test_chemistry.py` certifies H and He nuclei
  conservation at each level, charge conservation
  (n_e = HII + HeII/4 + HeIII/2), and physical direction — warm dense gas
  recombines and cools, cold gas forms H2 (9-species) and HD (12-species).

This is the foundation for wrapping chemistry and the photon transport solver
— they reuse the same primitives; see `docs/WRAPPING.md`.

- **MHD constrained transport** (`bridge.mhdct_step`) -- the CT solver
  `grid::SolveMHD_Li` (HydroMethod=MHD_Li, UseMHDCT).  `tests/test_mhdct.py`
  certifies the defining CT invariant on a Brio-Wu tube: the discrete
  divergence of the face-centered field stays at machine zero under evolution.
- **AMR control** (`bridge.flag_cells`, `bridge.cluster`) -- cell flagging
  via the real `grid::SetFlaggingField` dispatch (any CellFlaggingMethod), with
  the criteria certified: slope, baryon-mass/overdensity, second derivative,
  shear, must-refine region (geometric), and Jeans length (gravity)
  (`tests/test_amr.py`); plus Berger-Rigoutsos clustering
  of flagged cells into child grids (`ProtoSubgrid` +
  `IdentifyNewSubgridsBySignature`).  The remaining criteria (particle mass,
  cooling time, optical depth, resistive length, metallicity, metal mass)
  dispatch through the same path once their fields/globals are set.
- **Time integration** (`bridge.compute_timestep`, `problems.Problem(...,
  evolve=True)`) -- the CFL timestep (`grid::ComputeTimeStep`, certified =
  courant*dx/(|v|+c_s)) and Enzo's full AMR time integrator
  (`EvolveHierarchy`: per-grid timesteps -> boundary conditions -> the
  hydro/MHD/gravity solvers -> AMR sub-cycling -> flux correction ->
  projection -> hierarchy rebuild, looped to StopTime).  `tests/test_evolve.py`
  evolves the Toro-1 shock tube end-to-end and matches the exact Riemann
  solution (L1 ~ 3e-3).
- **AMR hierarchy operators** (`bridge.project_to_parent`,
  `bridge.interpolate_to_child`, `bridge.correct_refined_fluxes`) -- the
  inter-grid operators between refinement levels: restriction
  (`grid::ProjectSolutionToParentGrid`, fine->coarse averaging), prolongation
  (`grid::InterpolateFieldValues`, coarse->fine), and refluxing
  (`grid::CorrectForRefinedFluxes`).  `tests/test_hierarchy.py` certifies each
  by conservation/accuracy on a parent+child pair: restriction exact to 0,
  linear prolongation to 4e-16, and refluxing corrects exactly the coarse cells
  at the fine-grid boundary.

### Problem setup / initial conditions (all problem types)

`enzomodules.problems` drives Enzo's own `InitializeNew`, which reads a `.enzo`
parameter file and dispatches on `ProblemType` to **every** problem generator —
so one wrapper produces the initial conditions of *any* Enzo problem type and
exposes the resulting grid hierarchy:

```python
from enzomodules.problems import Problem
with Problem("run/Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo") as p:
    g = p.grid(0)                       # walk the hierarchy
    rho = g.field("Density")            # flat field data (incl. ghost zones)
    print(p.problem_type, g.rank, g.dims, g.field_names, g.num_particles)
```

`tests/test_problems.py` initializes a spread of problem types (Sod shock tube,
Sedov blast, Implosion, Kelvin-Helmholtz, Noh) through the single dispatch and
checks the generated fields — e.g. the Sod IC has the correct left/right
density states.  This is exactly what a rewrite needs to cross-check its own
initial-condition generators against the reference, problem by problem.

`tools/capture_problems.py` snapshots a compact **golden signature** of each
problem's initial conditions (per-field dims + min/max/mean/sum, per the
`fixtures/Problems/*.json`); `tests/test_problem_fixtures.py` re-initializes
and compares, catching any drift in Enzo's IC generators — and giving a rewrite
a reference to diff its own generators against.

> Self-contained (analytic) problems work directly; problem types that read
> external data (e.g. cosmological initial-condition files) need that data
> present, as in a normal Enzo run.

### Driving the time loop from Python (`problems.Session`)

`Problem(..., evolve=True)` hands the whole run to Enzo's monolithic
`EvolveHierarchy`.  A **`Session`** instead holds the *live* hierarchy and
exposes each orchestration step `EvolveLevel` runs, so a host language can own
the time loop — a Python (or Julia) re-implementation of `EvolveLevel` that
calls the certified legacy reference for each step:

```python
from enzomodules.problems import Session
with Session("run/Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo") as s:
    while s.time < s.stop_time:
        s.set_boundary(0)                 # SetBoundaryConditions (FAST_SIB path)
        dt = s.compute_dt(0); s.set_dt(0, dt)   # grid::ComputeTimeStep (CFL)
        s.solve_hydro(0)                  # grid::SolveHydroEquations
        s.advance_time(0)                 # SetTimeNextTimestep + cycle bump
    rho = s.grid(0).field("Density")      # read the evolving state
```

The steps exposed on the live hierarchy:

| `Session` method | legacy call |
| --- | --- |
| `set_boundary(level)` | `CreateSiblingList` + `SetBoundaryConditions` |
| `compute_dt(level)` / `set_dt` | `grid::ComputeTimeStep` / `SetTimeStep` |
| `solve_hydro(level)` | `grid::SolveHydroEquations` |
| `advance_time(level)` | `grid::SetTimeNextTimestep` + cycle bump |
| `gravity(level)` | `PrepareDensityField` (Poisson) + `ComputeAccelerations` + `CopyPotentialToBaryonField` |
| `update_particles(level)` | `UpdateParticlePositions` |
| `rebuild(level)` | `RebuildHierarchy` (flag + cluster + regrid) |
| `evolve_photons(level)` | `EvolvePhotons` |
| `create_fluxes` / `finalize_fluxes(level)` | `CreateFluxes` / `FinalizeFluxes` (per-level boundary-flux storage) |
| `clear_boundary_fluxes(level)` | `grid::ClearBoundaryFluxes` |
| `copy_baryon_to_old(level)` | `grid::CopyBaryonFieldToOldBaryonField` |
| `update_from_finer(level)` | `CreateSUBlingList` + `UpdateFromFinerGrids` (project + flux-correct) |
| `num_grids_on_level(level)` | grids per level (for the recursion) |
| `step(level)` / `run(level)` | the default root-level loop |

The headline guarantee (`tests/test_session.py`): the Python-driven loop
reproduces `EvolveHierarchy` **bit-for-bit** on the Toro-1 shock tube
(`Linf = 0`) and matches the exact Riemann solution — the equivalence that lets
`EvolveLevel` be re-implemented step-by-step on top of these certified bridges
and verified at every stage.  The gravity and particle steps run on the live
hierarchy; their physics is certified in the dedicated bridge tests
(`test_gravity.py`, `test_particles.py`).

`evolve_photons` is **fully emitting**: it builds the `SubgridMarker`
grid-ownership map, sizes the photon timestep, and sub-cycles the photon time so
the sources actually radiate — depositing photo-ionization / heating rates and
coupling to the chemistry.  `tests/test_session.py` drives the `PhotonTest`
Strömgren-sphere problem through the session and checks that the ionization
front propagates (the ionized fraction grows monotonically, the illuminated
volume expands) — the full `EvolvePhotons` orchestration, not just the
ray-tracer kernels (which are separately certified in `test_radiation.py`).

Two source modes:
- `evolve_photons(level)` — the radiation sources read from the parameter file
  (e.g. `PhotonTest`).
- `evolve_photons(level, stars=True)` — gather the grid's **star particles**
  (`StarParticleInitialize`) and convert the active ones into radiation sources
  (`StarParticleRadTransfer`), so star particles radiate.  A star only emits once
  Enzo's formation lifecycle has marked it an active radiation source (`type > 0`,
  alive, non-formation feedback flag) — which `star_particles` (below) drives.

**RT-stack hardening.**  Because this library keeps one set of Enzo globals for
the whole process, radiative-transfer state could leak between `Session`
instances.  `evolve_photons` is a clean no-op when `RadiativeTransfer` is off or
a level is empty, builds the `SubgridMarker` the photon walk dereferences (a
former segfault), sizes `dtPhoton` *without* calling `RadiativeTransferPrepare`
(whose `StarParticleRadTransfer` would wipe the parameter-file sources), and
`session_init` resets `RadiationPressure` (reset only on the RT-on path in stock
Enzo, so it otherwise leaked into a later pure-hydro session and made the PPM
solver dereference a missing gravity field).  `tests/test_session.py` exercises
the RT→non-RT transition to lock this down.

#### Writing EvolveLevel itself in Python

The last group of steps above (flux storage + `update_from_finer`) is the AMR
*conservation* machinery, which is exactly what makes `EvolveLevel` non-trivial.
With it, `Session.evolve_level()` **is** a Python re-implementation of Enzo's
recursive `EvolveLevel`, built entirely from the certified bridges — it clears
the boundary fluxes, sub-cycles the level to its parent's timestep, recurses
into the finer level between sub-steps, then projects the fine solution back and
corrects the coarse fluxes (`update_from_finer`) and regrids:

```python
with Session("run/Hydro/Hydro-2D/ImplosionAMR/ImplosionAMR.enzo") as s:
    s.run_amr()                       # initial regrid, then evolve_level(0) + regrid to StopTime
    print(s.grid(0).field("Density"))
```

`evolve_level` / `run_amr` take `gravity=True` (run the self-gravity chain each
step), `radiation=True` (emit + transport photons each step),
`star_sources=True` (radiate from star particles) and `cooling=True` (solve
radiative cooling + chemistry each step), so the same Python driver does
gravitating, radiation-hydrodynamic, cooling AMR runs.

`tests/test_session.py` runs this on the 4-level `ImplosionAMR` problem — the
recursion, sub-cycling, flux correction, projection and regridding all fire (the
hierarchy grows from `[1,1,1,1]` to hundreds of subgrids) — and the result
matches Enzo's own `EvolveHierarchy` on the root grid to `L1 ~ 4e-5`, with the
mean density (a conservation proxy) agreeing to six figures.  It also drives the
coupled radiation-hydrodynamics problem `PhotonTestAMR` (hydro + RT + AMR) with
`radiation=True`, ionizing the medium as the hydro and regridding run.  In other
words, the whole AMR time integrator — `enzo.C`, `EvolveHierarchy`, *and* the
recursive `EvolveLevel`, including self-gravity and radiative transfer — can now
be written in Python on top of EnzoModules, with each step verified against the
legacy reference.

#### Inline halo finder (FOF) + SUBFIND

`Session.find_halos(subfind=False, linking_length=0.0, min_size=0)` runs Enzo's
inline friends-of-friends halo finder — the same `FOF()` pipeline
`EvolveHierarchy` runs inline — on the session's live particle hierarchy, and
returns the catalogue *in memory* (a list of halo dicts: particle count, total /
virial mass, virial radius, centre of mass, mean velocity, velocity dispersion,
spin, angular momentum) sorted largest-first, instead of only writing it to
disk.  Pass `subfind=True` to also run SUBFIND; `subhalo_count()` then reports
the number of gravitationally self-bound subgroups.  Halo properties come from
the finder's own `get_particles`/`get_properties`, so the computation is the
certified one, not re-derived.

The finder *moves* the particles off the grids (`FOF_Initialize` /
`MoveParticlesFOF`) and the bridge calls `FOF_Finalize` to restore and
redistribute them, so `find_halos` is **repeatable and non-destructive**: the
total particle count is conserved and the session keeps evolving afterwards.  It
is a clean no-op (returns `[]`) on a particle-free problem, and `min_size` is
clamped to SUBFIND's neighbour requirement so a small value can't abort the run.
`tests/test_session.py` drives it on `GravityTest` (5000 particles): FOF groups
grow with the linking length, SUBFIND resolves subgroups in the big clump,
particle count is conserved across repeated calls, and the hydro still advances
afterwards.

#### Radiative cooling + non-equilibrium chemistry

`Session.solve_cooling(level)` runs `grid::MultiSpeciesHandler` — the
cooling/chemistry sub-step `EvolveLevel` runs after the hydro solve (it is a
*separate* step, not part of `solve_hydro`).  It dispatches to Grackle, the
coupled rate-and-cool solver, or `SolveRateEquations` + `SolveRadiativeCooling`
per the run's settings, and is a clean no-op when both `MultiSpecies` and
`RadiativeCooling` are off.  `tests/test_session.py` heats/ionizes the
`PhotonTest` gas with a few RT steps, then shows repeated `solve_cooling` (no
further heating) **monotonically lowers the total energy and evolves the
species** — real cooling/chemistry on the live hierarchy, not a callable no-op.

This was the top P0 gap in `docs/COVERAGE_AUDIT.md`; with it, the Python driver
can run a cooling, radiation-hydrodynamic, gravitating AMR step.  Wiring it also
hardened `session_init`: a problem whose initialization throws (e.g. a missing
cooling-rate data file) now raises a Python error instead of aborting the host
process.

#### Data output + checkpoint/restart

`Session.write_output(number)` writes the full state to disk
(`Group_WriteAllData`) — the same HDF5 dump Enzo produces: grid data, the
hierarchy, the external boundary and a parameter file — and returns the dump's
path.  `Session.from_output(path)` reloads it (`Group_ReadAllData`) into a fresh
session.  Together they are a true **checkpoint/restart**: `tests/test_session.py`
shows that running N steps, checkpointing, reloading and continuing M steps
reproduces an uninterrupted N+M-step run **bit-for-bit**, and that a multi-grid
AMR hierarchy round-trips (grid count + state preserved).

This was P0 gap #2.  It also closes a certification loop: dumps can now be diffed
against on-disk Enzo output.  Two robustness points fell out of the wiring:
Enzo's output-name buffer is only initialized when the basename matches a
recognized dump-name pattern (otherwise the name is built from uninitialized
stack memory), so the bridge passes the run's `DataDumpName` and writes directly
to the session's working directory with a deterministic `<DataDumpName><id>`
name; and `from_output` raises a Python error on a bad dump instead of aborting.

#### Star formation, feedback & activation (radiating stars)

`Session.star_particles(level)` runs the star-particle lifecycle EvolveLevel
drives: `StarParticleInitialize` (gather stars, set feedback flags) → per-grid
`StarParticleHandler` (star formation + feedback, `STARMAKE_METHOD`) →
`StarParticleFinalize` / `ActivateNewStar` (flip newly-formed / aged particles
into **active, radiating** stars and sync back to the grids).  No-op unless
`StarParticleCreation` or `StarParticleFeedback` is on.

This was P0 gap #3, and it closes the loop on radiating stars: previously a Pop
III star stayed *unborn* (`type < 0`), so `evolve_photons(stars=True)` was a
no-op.  With `star_particles` it **activates and radiates** —
`tests/test_session.py` drives `ProblemType 252` (a single Pop III star) with
`star_particles` + `evolve_photons(stars=True)` and shows the star ionizes the
surrounding gas into a Strömgren sphere whose illuminated volume grows
monotonically (a propagating I-front).  The one subtlety is temporal resolution:
a star radiates only during its (short, for Pop III) main-sequence lifetime, so
the timestep must be finer than that lifetime to catch it — with the default
light-crossing dt the star is born and supernovae within one step.  In a Python
driver: `run_amr(radiation=True, star_sources=True, star_formation=True)`.

#### More physics methods & modules

Rounding out the audit's P1/P2 items, the session also exposes:

- **All hydro/MHD methods** — `solve_hydro` dispatches the Runge-Kutta family
  (HD_RK / MHD_RK: two-step integration with Dedner wave speeds and `NColor`
  setup) in addition to PPM/Zeus.  `tests/test_session.py` runs the Brio-Wu MHD
  shock tube (MHD_RK) step-by-step, finite and mass-conserving.  (Constrained-
  transport MHD dispatches but needs extra field orchestration — see
  `docs/COVERAGE_AUDIT.md`.)
- **Active particles** — `active_particles(level)` (the sink / SmartStar /
  accretion framework: ActiveParticleInitialize → Handler → Finalize).
- **Cosmology** — `cosmology()` / `scale_factor` / `redshift`
  (CosmologyComputeExpansionFactor).
- **UV background** — `update_radiation_field(level)` (RadiationFieldUpdate).
- **Turbulence driving** — `random_forcing(level)`
  (ComputeRandomForcingNormalization + ComputeStochasticForcing).
- **Thermal conduction** (`conduct_heat`) and **shock finding** (`find_shocks`).
- **Diagnostics & hooks** — `domain_boundary_mass_flux(level)` and
  `problem_specific_routines(level)`.

Each is a clean no-op when its physics is off, so they compose freely into the
Python `evolve_level` / `run_amr` driver.  `docs/COVERAGE_AUDIT.md` tracks what
remains explicitly deferred (implicit FLD, CT-MHD orchestration, libyt).

- **Multi-grid (AMR) photon transport** (`bridge.raytrace_twogrid`) — a ray
  crossing two tiled grids (the legacy `SubgridMarker` -> `FindPhotonNewGrid`
  handoff) attenuates exactly as one grid of the combined length
  (`tests/test_radiation.py`).
- **Multi-dimensional PPM** (`bridge.ppm_hydro_step`) — the full
  `grid::SolveHydroEquations`, which runs the x/y/z Euler sweeps.
  `tests/test_ppm_hydro.py` certifies the y-sweep against the x-sweep by
  rotational symmetry (profiles identical to <1e-10) and matches the exact
  Riemann solution.

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
