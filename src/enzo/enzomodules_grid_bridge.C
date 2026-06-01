/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- grid-method solvers (ZEUS, ...)
/
/  PURPOSE:
/    Exercise grid:: solver *methods* (which read a constructed grid object,
/    not plain arrays) from external ctypes callers.  Uses the generic grid
/    fixture primitives (Grid_EnzoModulesFixture.C) to build a minimal grid,
/    fill it, call the method, and read the result back.  Globals are seeded
/    with Enzo's own SetDefaultGlobalValues so solver-specific parameters
/    (e.g. ZEUS artificial viscosity, tiny_number) get sane values; only the
/    few that matter for the test are overridden.
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

int SetDefaultGlobalValues(TopGridData &MetaData);

extern "C" {

/* One ZEUS hydro update of a 1D slice (operator-split finite-difference
 * solver, grid::ZeusSolver).  d = density, e = specific *internal* energy,
 * u = x-velocity, each length idim (= active + 2*nghost).  Updated in place.
 * Returns 0 on success. */
int enzomodules_zeus_sweep_1d(double *d, double *e, double *u,
                              int idim, int nghost, double dx, double dt,
                              double gamma)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);     /* sane defaults for every global */

  Gamma               = gamma;
  HydroMethod         = Zeus_Hydro;
  NumberOfGhostZones  = nghost;
  DualEnergyFormalism = 0;
  PressureFree        = 0;
  ComovingCoordinates = 0;
  UseDrivingField     = 0;

  grid g;
  int   dims[3]  = { idim, 1, 1 };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)(dx * idim), (FLOAT)1.0, (FLOAT)1.0 };
  int   ftypes[5] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };

  int size = g.EnzoModulesSetupGrid(1, dims, left, right, 5, ftypes, dt);

  int iDen = g.EnzoModulesFieldIndex(Density);
  int iTE  = g.EnzoModulesFieldIndex(TotalEnergy);
  int iV1  = g.EnzoModulesFieldIndex(Velocity1);
  int iV2  = g.EnzoModulesFieldIndex(Velocity2);
  int iV3  = g.EnzoModulesFieldIndex(Velocity3);
  std::vector<double> zero(size, 0.0);
  g.EnzoModulesSetField(iDen, d);
  g.EnzoModulesSetField(iTE,  e);
  g.EnzoModulesSetField(iV1,  u);
  g.EnzoModulesSetField(iV2,  zero.data());
  g.EnzoModulesSetField(iV3,  zero.data());

  float gam = (float)gamma;
  std::vector<float> dxa(idim, (float)dx), dya(1, (float)dx), dza(1, (float)dx);
  long_int GridGlobalStart[3] = { 0, 0, 0 };
  fluxes  *SubgridFluxes[1]   = { NULL };
  int      colnum[1]          = { 0 };

  int rc = g.ZeusSolver(&gam, /*igamfield*/ 0, /*nhy*/ 0,
                        dxa.data(), dya.data(), dza.data(),
                        /*gravity*/ 0, /*NumberOfSubgrids*/ 0, GridGlobalStart,
                        SubgridFluxes, /*NumberOfColours*/ 0, colnum,
                        /*bottom*/ 1, /*minsupecoef*/ 0.0);

  if (rc != FAIL) {
    g.EnzoModulesGetField(iDen, d);
    g.EnzoModulesGetField(iTE,  e);
    g.EnzoModulesGetField(iV1,  u);
  }
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
