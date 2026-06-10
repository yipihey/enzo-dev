# ADR-0006 — Multiphysics on the GPU: KA chemistry + M1 radiative transfer

Status: **accepted, partially implemented** (see "Status & staging" below).

## Context

EnzoNG's hydro (`PPMKernels`) and gravity (`PoissonKernels`) are single-source
KernelAbstractions kernels: one precision-generic `@kernel` runs on the CPU (the
f64 parity oracle, certified bit-tight against Enzo's Fortran) and on Metal
(f32). To "truly enable multiphysics on the GPU" we need the two remaining pillars
of Enzo's microphysics on the same footing:

1. **Non-equilibrium multi-species chemistry + radiative cooling** — the
   `calc_rates.F → cool1d_multi.F → solve_rate_cool.F` chain (Abel, Anninos,
   Zhang & Norman 1997).
2. **Radiative transfer** — coupling the radiation field back into the chemistry.

## Decision

### Chemistry → `lib/ChemistryKernels`

Port the network as one **size- and precision-generic** kernel, switchable at the
call site on `Val{NSP}` for the three Enzo `MultiSpecies` networks:

| `nspecies` | `MultiSpecies` | species |
|---|---|---|
| 6  | 1 | de HI HII HeI HeII HeIII |
| 9  | 2 | + HM H2I H2II |
| 12 | 3 | + DI DII HDI |

* Rate/cooling tables (`calc_rates.F`) are built **host-side once** (Abel+97
  collisional set, Cen-92 cooling, case-A/B recombination) and uploaded —
  identical to Enzo's table construction.
* The integrator (`solve_rate_cool.F`) is **per-cell**: one work-item owns one
  cell's entire sub-cycled stiff solve, using Anninos+97's semi-implicit
  backward-difference update `n_new = (creation·dt + n_old)/(1 + destruction·dt)`,
  with adaptive sub-stepping and a `make_consistent` renormalization (Grackle's
  step) so H and He nuclei are conserved exactly.
* This per-cell independence is *why* chemistry is the most natural GPU target in
  Enzo — the same shape as `pgas2d_dual!`'s row-sequential sweep.

### Radiative transfer → `lib/RadiationKernels` (M1, not ray tracing)

Enzo's production RT is **adaptive ray tracing** (Abel & Wandelt; HEALPix photon
packages traversing AMR with MPI transport). That is irregular pointer-chasing
across grids — it does **not** map to a flat KA kernel, and a faithful GPU port
would be a research project with limited efficiency.

Instead we implement the **two-moment M1** method (Aubert & Teyssier 2008;
Rosdahl+ 2013 / RAMSES-RT): carry `(Nγ, Fγ)` per photon group, close with the
Levermore Eddington tensor, transport with a dimensionally-split global
Lax-Friedrichs scheme, and couple per-cell into the chemistry network via
photo-ionization/heating (`RadiativeTransferIonization.C`). M1 is a stencil + a
per-cell update — the same GPU-efficient shape as the hydro/Poisson kernels — and
supports a reduced speed of light. Group count is a call-site `Val{NG}`.

The two libraries meet at `ChemistryKernels.PhotoRates`: `RadiationKernels`
fills `kphHI/kphHeI/kphHeII/kdissH2/photogamma`, which `solve_rate_cool!`
consumes as source terms.

## Consequences

* One source per physics, CPU + Metal, matching the existing kernel libraries
  (backend registry, `device_zeros`/`to_device`/`to_host`, Metal weakdep
  extension, `Adapt` rules for the struct-of-arrays containers).
* f32 on Metal is an accuracy compromise for stiff chemistry — the same one
  GPU-Grackle makes; the f64-CPU parity layer is where it gets policed.
* We diverge from Enzo's *algorithm* for RT (M1 vs ray tracing). The parity
  oracle for RT is therefore an M1 reference, not Enzo's ART — a deliberate
  choice recorded here.

## Status & staging

Implemented and structured for the parity harness (CPU correctness/conservation
suites included; they run headless on any CI):

- ✅ `ChemistryKernels`: rate/cooling tables, `cool1d_multi!` edot, switchable
  `solve_rate_cool!` for **6-, 9- AND 12-species** — the full collisional +
  H2/H- network plus the **deuterium/HD set (k50–k56, reactions 50–56)** with the
  DI/DII/HDI semi-implicit updates and Enzo's exact charge-conservation electron
  closure. `conserve` flag toggles `make_consistent` (default on for robustness;
  **off for bit-exact Enzo parity**, since native `solve_rate_cool` does not
  renormalize species).
- ✅ `RadiationKernels`: M1 closure, GLF multi-group transport, per-cell photo
  coupling.
- ✅ **Generic mass-action executor** (`ChemistryKernels.GenericNetwork` /
  `generic_step!`): an arbitrary-network sibling of the primordial solver, the
  lowering target for KROME networks (CPU+GPU).
- ✅ **`lib/Krome`** (Krome.jl): parses KROME `react_*` files, compiles their
  Fortran rate expressions to Julia, and lowers onto `GenericNetwork`. Upstream
  KROME is vendored as a **git submodule** (`extern/krome`, GPL-3.0, see the
  license note in `lib/Krome/README.md`).

Remaining (explicit stubs, next stages):

- ⏳ **Bit-tight Fortran-fixture certification.** The suites here are f64 sanity +
  conservation gates plus a gated parity scaffold
  (`test/test_parity_fortran.jl`). The golden-fixture layer diffs KA
  `solve_rate_cool!(conserve=false)` against Enzo's `solve_rate_cool` via
  `EnzoLib.session_solve_cooling` (same mechanism PPMKernels/PoissonKernels use);
  the remaining hookup is a one-zone field accessor on the bridge. Not runnable in
  the porting environment (no Julia runtime).
- ⏳ **KROME DLSODES-parity.** The generic executor uses the Anninos+97
  semi-implicit scheme, not KROME's DLSODES; density/dust/shielding-dependent
  rates are skipped (counted in `KromeNetwork.skipped`). Temperature-dependent
  networks (incl. all primordial variants) are covered.
- ⏳ **Cloudy/metal cooling** (`cool1d_cloudy.F`) and the **UV-background table**
  (`RadiationFieldCalculateRates.C`) — additive, not on the critical path.
