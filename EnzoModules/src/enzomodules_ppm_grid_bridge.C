/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- full multi-dimensional PPM hydro step
/
/  PURPOSE:
/    Exercise grid::SolveHydroEquations, which internally runs the x, y, and z
/    directional PPM Euler sweeps (Grid_xEulerSweep/yEulerSweep/zEulerSweep)
/    when HydroMethod == PPM_DirectEuler.  Uses the generic grid fixture
/    primitives (Grid_EnzoModulesFixture.C) to build a minimal grid, fill the
/    Density/TotalEnergy/Velocity1-3 fields, take one full multi-dim hydro
/    step in place over dtFixed, and read the fields back.  Globals are seeded
/    with Enzo's SetDefaultGlobalValues, then the PPM-relevant ones overridden.
/
/    Compiled with Enzo's headers and linked against the Enzo shared library.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

extern "C" void enzomodules_init_timer();

int SetDefaultGlobalValues(TopGridData &MetaData);

extern "C" {

/* One full multi-dimensional PPM (PPM_DirectEuler) hydro step of a grid of
 * the given rank/dims, via grid::SolveHydroEquations (which dispatches to the
 * x/y/z Euler sweeps).  The domain spans [0, dims[k]*dx] on each axis with a
 * uniform cell width dx.  Fields are Density, TotalEnergy (specific total
 * energy = internal + 0.5 v^2), Velocity1/2/3, each a flat row-major array of
 * length prod(dims) INCLUDING the ghost zones the caller laid out.  The grid
 * is advanced in place over dtFixed = dt and the updated fields are written
 * back to d/e/u/v/w.  Returns 0 on success. */
int enzomodules_ppm_hydro_step(int rank, int dims[3], double dx, double dt,
                               double gamma,
                               double *d, double *e,
                               double *u, double *v, double *w)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);     /* sane defaults for every global */
  enzomodules_init_timer();

  Gamma                      = gamma;
  HydroMethod                = PPM_DirectEuler;
  RiemannSolver              = TwoShock;
  RiemannSolverFallback      = 0;
  ReconstructionMethod       = PPM;
  ConservativeReconstruction = 0;
  PositiveReconstruction     = 0;
  DualEnergyFormalism        = 0;
  EOSType                    = 0;
  PressureFree               = 0;
  ComovingCoordinates        = 0;
  NumberOfGhostZones         = 3;
  UseHydro                   = 1;
  SelfGravity                = 0;
  UniformGravity             = 0;
  PointSourceGravity         = 0;
  UseMinimumPressureSupport  = FALSE;
  MultiSpecies               = 0;
  RadiativeTransfer          = 0;
  RadiativeTransferFLD       = 0;
  UseMHDCT                   = 0;
  CRModel                    = 0;
  ShockMethod                = 0;
  MaximumRefinementLevel     = 0;

  grid g;
  int   gdims[3] = { 1, 1, 1 };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  for (int dim = 0; dim < rank; dim++) {
    gdims[dim] = dims[dim];
    left[dim]  = (FLOAT)0.0;
    right[dim] = (FLOAT)(dx * dims[dim]);
  }
  int ftypes[5] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };

  int size = g.EnzoModulesSetupGrid(rank, gdims, left, right, 5, ftypes, dt);

  int iDen = g.EnzoModulesFieldIndex(Density);
  int iTE  = g.EnzoModulesFieldIndex(TotalEnergy);
  int iV1  = g.EnzoModulesFieldIndex(Velocity1);
  int iV2  = g.EnzoModulesFieldIndex(Velocity2);
  int iV3  = g.EnzoModulesFieldIndex(Velocity3);

  g.EnzoModulesSetField(iDen, d);
  g.EnzoModulesSetField(iTE,  e);
  g.EnzoModulesSetField(iV1,  u);
  g.EnzoModulesSetField(iV2,  v);
  g.EnzoModulesSetField(iV3,  w);

  int rc = g.SolveHydroEquations(0, 0, NULL, 0);

  if (rc != FAIL) {
    g.EnzoModulesGetField(iDen, d);
    g.EnzoModulesGetField(iTE,  e);
    g.EnzoModulesGetField(iV1,  u);
    g.EnzoModulesGetField(iV2,  v);
    g.EnzoModulesGetField(iV3,  w);
  }
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
