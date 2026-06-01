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
int MultigridSolver(float *RHS, float *Solution, int Rank, int TopDims[],
                    float &norm, float &mean, int start_depth,
                    float tolerance, int max_iter);

extern "C" {

/* Gravity/Poisson: solve the discrete Poisson equation L(phi) = rhs on a
 * uniform grid with Enzo's multigrid solver (the engine behind
 * grid::SolveForPotential).  dims are the field dimensions (use 2^k+1 per
 * axis for clean multigrid coarsening); rhs_in and solution_out are flat
 * length prod(dims).  Returns the converged residual diagnostics in
 * *out_norm / *out_mean.  Returns 0 on success. */
int enzomodules_poisson_solve(int rank, int dims[],
                              const double *rhs_in, double *solution_out,
                              double *out_norm, double *out_mean)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);

  int size = 1;
  for (int d = 0; d < rank; d++) size *= dims[d];

  float *rhs = new float[size];
  float *sol = new float[size];
  for (int i = 0; i < size; i++) { rhs[i] = (float)rhs_in[i]; sol[i] = 0.0; }

  float norm = 0.0, mean = 0.0;
  int tdims[3] = { 1, 1, 1 };
  for (int d = 0; d < rank; d++) tdims[d] = dims[d];

  int rc = MultigridSolver(rhs, sol, rank, tdims, norm, mean,
                           /*start_depth*/ 0, /*tolerance*/ 2.0e-6,
                           /*max_iter*/ 100);

  for (int i = 0; i < size; i++) solution_out[i] = (double)sol[i];
  *out_norm = (double)norm;
  *out_mean = (double)mean;
  delete[] rhs;
  delete[] sol;
  return (rc == FAIL) ? 1 : 0;
}

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

/* CIC-deposit particles onto a grid (the particle-mesh operation,
 * grid::DepositParticlePositions).  Builds a full grid with NumberOfParticles
 * particles at the given positions/masses, deposits to the gravitating mass
 * field, and returns that field (mass density) plus its cell count and cell
 * volume so callers can check CIC weights and mass conservation.
 * out_field must have capacity >= the deposit-field cell count; the actual
 * count is returned in *out_size.  Returns 0 on success. */
int enzomodules_cic_deposit(int rank, int dims[], double left[], double right[],
                            double *posx, double *posy, double *posz,
                            double *mass, int nparticles,
                            double *out_field, int out_capacity,
                            int *out_size, double *out_cellvol)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;
  SelfGravity         = 1;

  grid g;
  FLOAT le[3], re[3];
  for (int dim = 0; dim < 3; dim++) {
    le[dim] = (dim < rank) ? (FLOAT)left[dim]  : (FLOAT)0.0;
    re[dim] = (dim < rank) ? (FLOAT)right[dim] : (FLOAT)1.0;
  }
  int ftypes[1] = { Density };          /* one baryon field is enough */
  g.EnzoModulesSetupGrid(rank, dims, le, re, 1, ftypes, 0.0);

  g.EnzoModulesSetupParticles(nparticles, 0);
  g.EnzoModulesSetParticleMass(mass);
  g.EnzoModulesSetParticlePosition(0, posx);
  if (rank > 1) g.EnzoModulesSetParticlePosition(1, posy);
  if (rank > 2) g.EnzoModulesSetParticlePosition(2, posz);

  int size = g.EnzoModulesDepositParticles();
  if (size < 0 || size > out_capacity) { *out_size = size; return 1; }
  g.EnzoModulesGetDepositField(out_field);
  *out_size    = size;
  *out_cellvol = g.EnzoModulesDepositCellVolume();
  return 0;
}

/* Radiation transport: build a full grid carrying the radiative-transfer rate
 * fields (kphHI photo-ionization, PhotoGamma photo-heating), run Enzo's field
 * identification, and round-trip the kphHI field.  Proves the full-grid
 * fixture supports radiation-transport data structures and Enzo's RT field
 * machinery (the foundation for wrapping the photon-transport solver).
 * kph_in/kph_out are length idim; returns 0 on success and writes the located
 * field indices. */
int enzomodules_rt_identify(int idim, int nghost, double dx,
                            const double *kph_in, double *kph_out,
                            int *kphHINum_out, int *gammaNum_out)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  RadiativeTransfer             = 1;
  RadiativeTransferHydrogenOnly = 1;   /* only kphHI + PhotoGamma needed */
  MultiSpecies                  = 1;
  NumberOfGhostZones            = nghost;

  grid g;
  int   dims[3]  = { idim, 1, 1 };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)(dx * idim), (FLOAT)1.0, (FLOAT)1.0 };
  int   ftypes[3] = { Density, kphHI, PhotoGamma };
  g.EnzoModulesSetupGrid(1, dims, left, right, 3, ftypes, 0.0);

  int ikph = g.EnzoModulesFieldIndex(kphHI);
  g.EnzoModulesSetField(ikph, kph_in);

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int rc = g.IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                             kphHeIINum, kdissH2INum, kphHMNum,
                                             kdissH2IINum);
  *kphHINum_out = kphHINum;
  *gammaNum_out = gammaNum;
  g.EnzoModulesGetField(kphHINum, kph_out);
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
