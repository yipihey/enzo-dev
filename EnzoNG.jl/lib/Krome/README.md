# Krome.jl

A Julia front-end for the **KROME** astrochemistry package
([Grassi et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.439.2386G);
<https://bitbucket.org/tgrassi/krome>). It parses KROME `react_*` network files,
compiles their Fortran rate expressions to Julia, and lowers the network onto
[`ChemistryKernels.GenericNetwork`](../ChemistryKernels) so **any KROME network
runs on the CPU and the GPU** through KernelAbstractions — the same single-source
contract as EnzoNG's hand-written primordial network.

```
react_* file ──parse_krome──▶ KromeNetwork ──to_generic──▶ GenericNetwork ──generic_step!──▶ CPU/GPU
```

```julia
using Krome, ChemistryKernels
net = parse_krome(krome_network_path("react_primordial"))
gen = to_generic(Float64, net)                 # tabulated, backend-agnostic
be  = ChemistryKernels.backend(:cpu)           # or :metal
n   = ChemistryKernels.device_zeros(be, Float64, (gen.nspec, ncells))  # number densities, cm^-3
C, D = ChemistryKernels.allocate_workspace(be, Float64, gen.nspec, ncells)
ChemistryKernels.generic_step!(n, Tgas, gen, C, D; dt = dt_seconds)
```

## License — IMPORTANT

Upstream **KROME is GPL-3.0**. It is included here as a **git submodule**
(`extern/krome`), i.e. a *reference* to the upstream repository under its own
license — it is **not** copied into the EnzoNG tree, and EnzoNG's own code is not
relicensed. Krome.jl reads KROME's reaction databases (`networks/react_*`) as
data. If you redistribute or link KROME-derived Fortran, review the GPL-3.0
obligations and **cite Grassi et al. (2014)**.

Initialise the submodule before using the network-file features:

```
git submodule update --init EnzoNG.jl/lib/Krome/extern/krome
```

## Scope / limitations (current)

- **Temperature-dependent rates only.** Rate expressions in the supported
  `Tgas/Te/lnTe/invTe/…` variables compile; reactions whose rates reference local
  densities, dust, shielding, or user variables are **skipped** and counted in
  `KromeNetwork.skipped` (they need KROME's runtime context). The primordial
  networks are fully covered.
- **cgs units.** The generic executor uses number densities (cm⁻³) and seconds —
  distinct from `solve_rate_cool!`, which uses Enzo code units.
- **Mass-action semi-implicit integrator.** The executor uses the same
  Anninos+97 backward-difference scheme as the certified primordial solver, not
  KROME's DLSODES. For stiff networks needing tight tolerances, treat results as
  a fast approximate solver pending the DLSODES-parity layer.
