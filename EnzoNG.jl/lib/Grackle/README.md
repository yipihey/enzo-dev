# Grackle.jl

A Julia front-end for the **Grackle** chemistry & cooling library
([grackle-project/grackle](https://github.com/grackle-project/grackle), **Enzo
Public License**). Grackle descends from Enzo's `solve_rate_cool.F`/`cool1d_multi.F`
— the same primordial network already in [`ChemistryKernels`](../ChemistryKernels)
— so the new content wired here is Grackle's **Cloudy metal-line cooling**, ported
to KernelAbstractions as the GPU hot path.

```julia
using Grackle, ChemistryKernels
tbl = load_cloudy_table("…/solar_2008_3D_metals.h5"; redshift=0.0)   # → CloudyTable
# metal cooling on CPU or GPU, standalone:
ChemistryKernels.cloudy_metal_cooling!(edot_met, rhoH, logtem, Z, tbl; dom, idim, i1, i2)
# or integrated into the network energy update:
ChemistryKernels.solve_rate_cool!(ge, dens, sp, rt, u; nspecies, gamma, temperature_units,
    dt, idim, i1, i2, cloudy=tbl, metallicity=Z, fh=0.76)
```

## License

Grackle is under the **Enzo Public License** (BSD-style) — same lineage as EnzoNG,
no copyleft concern (unlike the GPL-3.0 KROME). Vendored as a git submodule under
`extern/grackle`. Please **cite the Grackle method paper** (Smith et al. 2017).
Initialise: `git submodule update --init EnzoNG.jl/lib/Grackle/extern/grackle`.

## Layout

- `cloudy_io.jl` — load Grackle Cloudy HDF5 tables (`/Parameter1`=log n_H,
  `/Temperature`, `/Cooling`,`/Heating`) into `ChemistryKernels.CloudyTable`.
  Rank-3 (n_H, z, T) tables are sliced at a redshift plane (full 3-D z
  interpolation is a planned extension).
- `ChemistryKernels/src/cloudy.jl` — the KA interpolation kernel
  (`cool1d_cloudy_g.F` port), CPU + Metal.
- `deps/build_grackle.sh` — build `libgrackle.so` for the binding/oracle
  (`grackle_version`, and the path to full `solve_chemistry` parity).

## Status

- ✅ Cloudy metal-cooling KA kernel — bit-exact vs the table at grid nodes,
  physical shape (peaks ~10⁵ K), linear in metallicity; integrated into
  `solve_rate_cool!`.
- ✅ `libgrackle` binding (version 3.4.2-dev verified).
- ⏳ Full end-to-end `solve_chemistry` parity (drive the C `chemistry_data` API
  from Julia) and rank-3 redshift interpolation — next steps.
