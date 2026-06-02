/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- chemistry & cooling (non-Grackle MultiSpecies)
/
/  PURPOSE:
/    Exercise Enzo's primordial chemistry + radiative cooling network
/    (grid::SolveRateAndCoolEquations -> solve_rate_cool) in isolation, for any
/    MultiSpecies level:
/        1 ->  6 species (e-, HI, HII, HeI, HeII, HeIII)
/        2 ->  9 species (+ HM, H2I, H2II)
/        3 -> 12 species (+ DI, DII, HDI)
/    Rate coefficients are computed analytically by InitializeRateData
/    (calc_rates), so no external data files are needed.
/
/    Builds a uniform block of gas with the requested species set, advances one
/    chemistry+cooling step, and returns the updated species and energy.
/    Compiled with Enzo headers, linked against libenzo.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

int SetDefaultGlobalValues(TopGridData &MetaData);
int InitializeRateData(FLOAT Time);

/* Canonical species ordering by MultiSpecies level. */
static int species_field_types(int multispecies, int *types)
{
  int k = 0;
  types[k++] = ElectronDensity;
  types[k++] = HIDensity;
  types[k++] = HIIDensity;
  types[k++] = HeIDensity;
  types[k++] = HeIIDensity;
  types[k++] = HeIIIDensity;
  if (multispecies > 1) {
    types[k++] = HMDensity;
    types[k++] = H2IDensity;
    types[k++] = H2IIDensity;
  }
  if (multispecies > 2) {
    types[k++] = DIDensity;
    types[k++] = DIIDensity;
    types[k++] = HDIDensity;
  }
  return k;   /* number of species */
}

extern "C" {

int enzomodules_num_species(int multispecies)
{ return (multispecies > 2) ? 12 : (multispecies > 1) ? 9 : 6; }

/* Advance one chemistry+cooling step on `n` cells at MultiSpecies level
 * `multispecies` (1|2|3).  rho/e_tot are length n; species_flat is
 * [nspecies][n] (row-major) in the canonical order above (nspecies from
 * enzomodules_num_species).  All densities are code-unit mass densities;
 * units are physical cgs scalings.  species_flat and e_tot are updated in
 * place.  Returns 0 on success. */
int enzomodules_chemistry_step(
    int multispecies, int n, double dt,
    double density_units, double length_units, double time_units,
    double *rho, double *e_tot, double *species_flat)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);

  MultiSpecies        = multispecies;
  RadiativeCooling    = 1;
  RadiationFieldType  = 0;
  ComovingCoordinates = 0;
  DualEnergyFormalism = 0;
  HydroMethod         = 0;
  Gamma               = 5.0 / 3.0;

  GlobalDensityUnits = density_units;
  GlobalLengthUnits  = length_units;
  GlobalTimeUnits    = time_units;
  GlobalMassUnits    = density_units * length_units * length_units * length_units;

  if (InitializeRateData(0.0) == FAIL) return 2;

  int sptypes[12];
  int nspecies = species_field_types(multispecies, sptypes);

  /* Field list: Density, TotalEnergy, Velocity1-3, then the species. */
  std::vector<int> ftypes;
  ftypes.push_back(Density);
  ftypes.push_back(TotalEnergy);
  ftypes.push_back(Velocity1);
  ftypes.push_back(Velocity2);
  ftypes.push_back(Velocity3);
  for (int s = 0; s < nspecies; s++) ftypes.push_back(sptypes[s]);

  grid g;
  int   dims[3]  = { n, 1, 1 };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  g.EnzoModulesSetupGrid(1, dims, left, right, (int)ftypes.size(), ftypes.data(), dt);

  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density),     rho);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(TotalEnergy), e_tot);
  for (int s = 0; s < nspecies; s++)
    g.EnzoModulesSetField(g.EnzoModulesFieldIndex(sptypes[s]), species_flat + (size_t)s * n);

  if (g.SolveRateAndCoolEquations(0) == FAIL) return 1;

  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(TotalEnergy), e_tot);
  for (int s = 0; s < nspecies; s++)
    g.EnzoModulesGetField(g.EnzoModulesFieldIndex(sptypes[s]), species_flat + (size_t)s * n);
  return 0;
}

} /* extern "C" */
