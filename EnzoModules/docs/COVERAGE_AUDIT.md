# EnzoModules coverage audit — toward a complete Python-callable Enzo

*What is exposed today, and what is still needed to drive a full Enzo
simulation (setup → a complete physics step → output/restart → analysis)
entirely from Python.*

The reference for "a complete step" is `EvolveLevel.C` (the recursive AMR time
integrator) plus the setup (`InitializeNew`) and teardown (`WriteAllData`) that
bracket it. Each item below is marked ✅ exposed as a live-hierarchy `Session`
step, 🟡 partial / kernel-level only, or ❌ missing.

---

## 1. What is already covered

### Session step surface (`problems.Session`, live hierarchy)
- ✅ **Setup** — `Session(paramfile)` → `InitializeNew`: every problem type, ICs
  from the parameter file (grid + particle problems).
- ✅ **Boundary conditions** — `set_boundary` (`SetBoundaryConditions`,
  interpolate from parent + external).
- ✅ **Timestep** — `compute_dt`, `set_dt` (CFL min over grids).
- ✅ **Hydro solve** — `solve_hydro` (`Grid_SolveHydroEquations`: PPM + Zeus).
- ✅ **Time / bookkeeping** — `advance_time`, `time`, `stop_time`, `cycle`.
- ✅ **Particles** — `update_particles` (drift/push on the live hierarchy).
- ✅ **Self-gravity** — `gravity` (`PrepareDensityField` + Poisson + accelerations).
- ✅ **AMR conservation machinery** — `rebuild` (`RebuildHierarchy`),
  `create_fluxes`, `update_from_finer` (flux correction + projection),
  `finalize_fluxes`, `clear_boundary_fluxes`, `num_grids_on_level`.
- ✅ **Radiative transfer (photon packages / Moray)** — `evolve_photons`
  (`EvolvePhotons`), parameter-file sources **and** star-particle sources
  (`evolve_photons(stars=True)`).
- ✅ **Inline halo finding** — `find_halos` (FOF) + SUBFIND, with the in-memory
  catalogue and `subhalo_count`.
- ✅ **Radiative cooling + non-equilibrium chemistry** — `solve_cooling`
  (`grid::MultiSpeciesHandler`: Grackle / coupled rate-and-cool / rate +
  cooling). *(was P0 gap #1; done)*
- ✅ **Drivers** — `evolve_level` (recursive Python EvolveLevel), `run_amr`,
  `step`, `run`.
- ✅ **Inspection** — field + particle getters per grid.

### Kernel-level (certified, but not yet wired as a live-hierarchy step)
- 🟡 `twoshock`, `ppm_sweep_1d` — PPM building blocks.
- 🟡 `hydro_rk_line`, `mhd_rk_line` — MUSCL/RK Riemann line solvers.
- 🟡 `mhdct_step` — constrained-transport MHD.
- 🟡 `chemistry_step` (`SolveRateEquations`) — non-equilibrium chemistry kernel.
- 🟡 `raytrace_*`, `hi_cross_section` — RT ray-tracer kernels.
- 🟡 AMR flagging (`flag_cells`, `flag_jeans`, `flag_region`), `interpolate_to_child`,
  `project_to_parent`, `correct_refined_fluxes`.

---

## 2. Gaps, by priority

### P0 — required for a physically complete step

| Gap | Enzo entry point | Why it matters | Suggested API |
|-----|------------------|----------------|---------------|
| ✅ ~~**Radiative cooling + non-equilibrium chemistry**~~ **— DONE** | `Grid::MultiSpeciesHandler` / `SolveRadiativeCooling` / `GrackleWrapper`, called from `EvolveLevel.C:671` — a **separate** step, *not* inside `solve_hydro` | The Session can do adiabatic + RT hydro but **cannot cool or evolve species** standalone. This was the single biggest missing physics step. | ✅ `session.solve_cooling(level)` (+ `cooling=True` in `evolve_level`/`run_amr`) |
| **Star formation + feedback** | `StarParticleInitialize` → star makers (`STARMAKE_METHOD`) → `Grid::StarParticleHandler` (feedback) → `StarParticleFinalize` / `ActivateNewStar` | Only `StarParticleInitialize` is reached today (inside `evolve_photons(stars=True)`). The full lifecycle — formation, feedback deposition, activation — is unexposed; this is why a Pop III star never *activates* to radiate. | `session.star_particles(level)` (init + handler + finalize) |
| **Data output / checkpoint + restart** | `WriteAllData`, `Group_WriteAllData`; `ReadAllData` | The Session cannot persist or reload state — no checkpoints, no handoff to yt on disk, no resume. | `session.write_output(name)`, `Session.from_output(dir)` |

### P1 — broaden physics coverage

| Gap | Enzo entry point | Notes | Suggested API |
|-----|------------------|-------|---------------|
| **Full hydro/MHD methods on the live hierarchy** | `Grid::*` RK solvers (`HD_RK`, `MHD_RK`), CT-MHD (`MHD_Li`), `ComputeDednerWaveSpeeds` | `solve_hydro` covers PPM/Zeus; the RK + CT-MHD live-hierarchy solves exist only as kernel line-solvers. | extend `solve_hydro` to dispatch all `HydroMethod`s |
| **Active particles** | `ActiveParticleInitialize` / `ActiveParticleFinalize` | Modern sink / SmartStar / accretion framework (`EvolveLevel` calls these alongside the star-particle path). | `session.active_particles(level)` |
| **Cosmology expansion** | `CosmologyComputeExpansionFactor`, comoving source terms | Implicit inside the solvers; no explicit expansion step / scale-factor accessor for a Python loop. | `session.cosmology_a()`, expansion in `advance_time` |
| **FLD radiation + UV background** | `RadiativeTransferCallFLD`, `RadiationFieldUpdate` | `evolve_photons` is the ray-tracing (Moray) RT only; implicit FLD and the homogeneous UV background are separate. | `session.solve_fld(level)`, `session.update_radiation_field()` |

### P2 — specialized physics & infrastructure

- ❌ **Stochastic forcing / turbulence driving** — `ComputeStochasticForcing`,
  `ComputeRandomForcingNormalization`.
- ❌ **Thermal conduction**, **cosmic rays**, **shock finding**, **MHD div-B
  cleaning** as explicit steps.
- ❌ **Conservation / diagnostics** — `CheckEnergyConservation`,
  `ComputeDomainBoundaryMassFlux`.
- ❌ **Problem-specific hooks** — `CallProblemSpecificRoutines`, event hooks.
- ❌ **Inline analysis** — `CallPython` / libyt streaming (halo finding is now
  done via `find_halos`).
- 🟡 **Gravity completeness** — `gravity` runs the chain; potential/acceleration
  boundary (`SetAccelerationBoundary`), external/static gravity and the APM
  subgrid solver should be checked for full multi-level coverage.

---

## 3. Certification gaps (orthogonal to feature coverage)

Exposing a step is necessary but not sufficient — EnzoModules' standard is a
test that pins the bridge against the legacy reference. Current state:

- ✅ **PPM** — bitwise vs `EvolveHierarchy`; exact-Riemann L1.
- ✅ **AMR integrator, gravity chain, RT (PhotonTest I-front), halo finder,
  cooling/chemistry** — behavioural tests (`solve_cooling` monotonic-cooling +
  chemistry test; graceful init-failure test).
- 🟡 **Zeus, RK hydro/MHD, CT-MHD** — kernel tests only; no live-hierarchy
  `solve_hydro` certification.
- ❌ **Star formation, cosmology, output/restart** — untested at the Session
  level (mostly because unexposed).

---

## 4. Recommended sequence

1. ✅ **`solve_cooling`** (P0) — **done**; unlocks realistic ISM/cosmology runs,
   pairs with the existing `chemistry_step` kernel and RT coupling.
2. **`write_output` + `from_output`** (P0) — checkpoint/restart; also lets tests
   diff against on-disk Enzo dumps for stronger certification.
3. **`star_particles`** (P0) — completes the star lifecycle and makes
   `evolve_photons(stars=True)` actually radiate.
4. **Extend `solve_hydro` to all `HydroMethod`s** (P1) — turns the kernel-level
   RK/CT-MHD work into live-hierarchy steps, with certification.
5. **`active_particles`, FLD/UV background, cosmology accessors** (P1).
6. **P2 specialized physics** as needed by target science problems.

With P0 complete, a Python driver could run a *cooling, star-forming,
radiating* (radiation-hydro + gravity + AMR) simulation end to end and check it
against a native Enzo run — the practical definition of "complete coverage."
