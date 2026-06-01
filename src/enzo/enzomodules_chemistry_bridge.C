/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- chemistry & cooling (non-Grackle MultiSpecies)
/
/  PURPOSE:
/    Exercise Enzo's primordial chemistry + radiative cooling network
/    (grid::SolveRateAndCoolEquations -> solve_rate_cool) in isolation.  The
/    rate coefficients are computed analytically by InitializeRateData
/    (calc_rates), so no external data files are needed.
/
/    Builds a uniform cell of gas with the 6-species set (e-, HI, HII, HeI,
/    HeII, HeIII), advances one chemistry+cooling step, and returns the updated
/    species and energy.  Compiled with Enzo headers, linked against libenzo.
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

extern "C" {

/* Advance one primordial-chemistry + cooling step on a uniform block of `n`
 * cells.  All species are mass densities in code units (HI/HII in H-mass,
 * He* in He-mass); `e_tot` is total specific energy.  Arrays length n are
 * updated in place.  Units are physical cgs scalings.  Returns 0 on success. */
int enzomodules_chemistry_step(
    int n, double dt,
    double density_units, double length_units, double time_units,
    double *rho, double *e_tot,
    double *de, double *HI, double *HII,
    double *HeI, double *HeII, double *HeIII)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);

  MultiSpecies        = 1;
  RadiativeCooling    = 1;
  RadiationFieldType  = 0;     /* no UV background */
  ComovingCoordinates = 0;
  DualEnergyFormalism = 0;
  HydroMethod         = 0;
  Gamma               = 5.0 / 3.0;

  GlobalDensityUnits = density_units;
  GlobalLengthUnits  = length_units;
  GlobalTimeUnits    = time_units;
  GlobalMassUnits    = density_units * length_units * length_units * length_units;

  if (InitializeRateData(0.0) == FAIL) return 2;

  grid g;
  int   dims[3]  = { n, 1, 1 };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  int   ftypes[11] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3,
                       ElectronDensity, HIDensity, HIIDensity,
                       HeIDensity, HeIIDensity, HeIIIDensity };
  g.EnzoModulesSetupGrid(1, dims, left, right, 11, ftypes, dt);

  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density),         rho);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(TotalEnergy),     e_tot);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(ElectronDensity), de);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HIDensity),       HI);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HIIDensity),      HII);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIDensity),      HeI);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIIDensity),     HeII);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIIIDensity),    HeIII);

  if (g.SolveRateAndCoolEquations(0) == FAIL) return 1;

  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(TotalEnergy),     e_tot);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(ElectronDensity), de);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(HIDensity),       HI);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(HIIDensity),      HII);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(HeIDensity),      HeI);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(HeIIDensity),     HeII);
  g.EnzoModulesGetField(g.EnzoModulesFieldIndex(HeIIIDensity),    HeIII);
  return 0;
}

} /* extern "C" */
