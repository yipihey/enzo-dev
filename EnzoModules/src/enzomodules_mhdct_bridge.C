/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- constrained-transport (CT) MHD step
/
/  PURPOSE:
/    Exercise Enzo's constrained-transport MHD solver, grid::SolveMHD_Li
/    (HydroMethod = MHD_Li, UseMHDCT = 1).  SolveHydroEquations dispatches to
/    SolveMHD_Li when UseMHDCT==1 && HydroMethod==MHD_Li, then computes the CT
/    electric field, takes the curl to update the FACE-CENTERED magnetic field
/    (MagneticField[]) and re-centers it.  This bridge builds a minimal grid
/    via the generic fixture, lays down both the cell-centered B (Bfield1/2/3)
/    and the staggered face-centered MagneticField[], measures the discrete
/    divergence of B before and after one CT step, and reads the fields back.
/    The defining CT property is that max|divB| stays at machine precision.
/
/    Mirrors enzomodules_ppm_grid_bridge.C: SetDefaultGlobalValues, override
/    the MHDCT-relevant globals, enzomodules_init_timer(), EnzoModulesSetupGrid
/    (which allocates the CT fields because UseMHDCT is set first), fill fields,
/    call g.SolveHydroEquations(0,0,NULL,0), read back.
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

/* One full constrained-transport MHD step (grid::SolveMHD_Li via
 * grid::SolveHydroEquations) of a 3D grid of the given dims (including ghost
 * zones).  The domain spans [0, dims[k]*dx] on each axis with uniform width dx.
 * d/e/u/v/w/bx/by/bz are flat row-major arrays of length prod(dims): density,
 * specific total energy, Velocity1/2/3, and the (cell-centered == face value
 * for a 1D-varying field) magnetic field components.  The face-centered
 * MagneticField[] is initialized from bx/by/bz directly (for a field that
 * varies only along x and is uniform transverse, the x-faces equal the cell
 * values, so this is divergence-free to machine precision).  Returns the max
 * |divB| over active cells before and after the step via the out-pointers, and
 * writes the updated fields back.  Returns 0 on success, 1 on solver failure. */
int enzomodules_mhdct_step(int dims[3], double dx, double dt, double gamma,
                           int nsteps,
                           double *d, double *e,
                           double *u, double *v, double *w,
                           double *bx, double *by, double *bz,
                           double *maxdivb_before, double *maxdivb_after)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);     /* sane defaults for every global */
  enzomodules_init_timer();

  Gamma                      = gamma;
  HydroMethod                = MHD_Li;
  UseMHDCT                   = 1;
  UseMHD                     = 1;        /* IdentifyPhysicalQuantities needs B */
  MaxVelocityIndex           = 3;        /* MHD always carries 3 velocities   */

  /* CT / MHD_Li solver parameters (see ReadParameterFile defaults for MHD_Li
     and the BrioWu-MHD-1D-MHDCT example). */
  RiemannSolver              = HLL;      /* robust for shock tubes (Brio-Wu)  */
  RiemannSolverFallback      = 0;
  ReconstructionMethod       = PLM;      /* MHD_Li default                    */
  MHD_CT_Method              = CT_Athena_LF;   /* 2: the default CT method    */
  MHDCTSlopeLimiter          = 1;
  MHDCTDualEnergyMethod      = 0;        /* DualEnergyFormalism == 0          */
  MHDCTPowellSource          = 0;
  MHDCTUseSpecificEnergy     = 1;        /* TE stored specific (we set it so) */
  MHD_ProjectE               = TRUE;     /* curl E -> update face B           */
  MHD_ProjectB               = FALSE;
  Theta_Limiter              = 1.0;
  EquationOfState            = 0;        /* adiabatic                          */
  IsothermalSoundSpeed       = 1.0;

  DualEnergyFormalism        = 0;
  ConservativeReconstruction = 0;
  PositiveReconstruction     = 0;
  PressureFree               = 0;
  EOSType                    = 0;
  ComovingCoordinates        = 0;
  NumberOfGhostZones         = 5;        /* MHD_Li / Brio-Wu use wide stencil */
  UseHydro                   = 1;
  SelfGravity                = 0;
  UniformGravity             = 0;
  PointSourceGravity         = 0;
  UseMinimumPressureSupport  = FALSE;
  MultiSpecies               = 0;
  RadiativeTransfer          = 0;
  RadiativeTransferFLD       = 0;
  CRModel                    = 0;
  ShockMethod                = 0;
  MaximumRefinementLevel     = 0;

  grid g;
  int   gdims[3] = { dims[0], dims[1], dims[2] };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)(dx*dims[0]), (FLOAT)(dx*dims[1]), (FLOAT)(dx*dims[2]) };

  int ftypes[8] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3,
                    Bfield1, Bfield2, Bfield3 };

  /* UseMHDCT is set above, so EnzoModulesSetupGrid -> PrepareGrid ->
     MHD_SetupDims and AllocateGrids allocate the face-centered MagneticField,
     ElectricField, AvgElectricField. */
  g.EnzoModulesSetupGrid(3, gdims, left, right, 8, ftypes, dt);

  int iDen = g.EnzoModulesFieldIndex(Density);
  int iTE  = g.EnzoModulesFieldIndex(TotalEnergy);
  int iV1  = g.EnzoModulesFieldIndex(Velocity1);
  int iV2  = g.EnzoModulesFieldIndex(Velocity2);
  int iV3  = g.EnzoModulesFieldIndex(Velocity3);
  int iB1  = g.EnzoModulesFieldIndex(Bfield1);
  int iB2  = g.EnzoModulesFieldIndex(Bfield2);
  int iB3  = g.EnzoModulesFieldIndex(Bfield3);

  g.EnzoModulesSetField(iDen, d);
  g.EnzoModulesSetField(iTE,  e);
  g.EnzoModulesSetField(iV1,  u);
  g.EnzoModulesSetField(iV2,  v);
  g.EnzoModulesSetField(iV3,  w);
  g.EnzoModulesSetField(iB1,  bx);
  g.EnzoModulesSetField(iB2,  by);
  g.EnzoModulesSetField(iB3,  bz);

  /* Fill the face-centered MagneticField[].  The faces normal to axis `dim`
     are an extra slab along that axis (MagneticDims = GridDimension + e_dim).
     For a Brio-Wu state the field varies only along x and is uniform in y,z:
     - Bx (x-faces): equals the cell value (Bx is constant), so we replicate
       the cell value into every x-face including the extra right face.
     - By, Bz (y-/z-faces): the cell value is exactly the face-centered value
       (no extra zone along the varying direction), so a direct copy of the
       per-cell array into the staggered slab matches at every shared face. */
  int gd0 = g.EnzoModulesMagneticSize(0); /* (nx+1)*ny*nz */
  int gd1 = g.EnzoModulesMagneticSize(1); /* nx*(ny+1)*nz */
  int gd2 = g.EnzoModulesMagneticSize(2); /* nx*ny*(nz+1) */
  int nx = dims[0], ny = dims[1], nz = dims[2];

  double *fx = new double[gd0];
  double *fy = new double[gd1];
  double *fz = new double[gd2];

  /* Bx on x-faces: MagneticDims[0] = {nx+1, ny, nz}. Cell (i,j,k) maps to two
     x-faces i and i+1; for a divergence-free 1D state Bx is x-constant, so we
     set face i from cell i and the extra right face nx from cell nx-1. */
  for (int k = 0; k < nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i <= nx; i++) {
        int ci = (i < nx) ? i : (nx - 1);
        int cell = ci + nx * (j + ny * k);
        int face = i + (nx + 1) * (j + ny * k);
        fx[face] = bx[cell];
      }

  /* By on y-faces: MagneticDims[1] = {nx, ny+1, nz}. By varies only along x;
     for each (i,k) every y-face equals that column's cell value. */
  for (int k = 0; k < nz; k++)
    for (int j = 0; j <= ny; j++)
      for (int i = 0; i < nx; i++) {
        int cj = (j < ny) ? j : (ny - 1);
        int cell = i + nx * (cj + ny * k);
        int face = i + nx * (j + (ny + 1) * k);
        fy[face] = by[cell];
      }

  /* Bz on z-faces: MagneticDims[2] = {nx, ny, nz+1}. */
  for (int k = 0; k <= nz; k++)
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
        int ck = (k < nz) ? k : (nz - 1);
        int cell = i + nx * (j + ny * ck);
        int face = i + nx * (j + ny * k);
        fz[face] = bz[cell];
      }

  g.EnzoModulesSetMagneticField(0, fx);
  g.EnzoModulesSetMagneticField(1, fy);
  g.EnzoModulesSetMagneticField(2, fz);

  delete [] fx;
  delete [] fy;
  delete [] fz;

  /* Initial divergence (over the persistent face-centered field). */
  *maxdivb_before = g.EnzoModulesMaxDivB();

  /* Take nsteps CT steps on the SAME grid so the face-centered MagneticField[]
     -- the true CT state -- is carried over natively between steps (only the
     cell-centered B is round-tripped through the C ABI, which is lossy, so we
     must iterate internally to test multi-step div-B preservation). */
  int rc = SUCCESS;
  double maxafter = *maxdivb_before;
  int ns = (nsteps < 1) ? 1 : nsteps;
  for (int step = 0; step < ns && rc != FAIL; step++) {
    /* Refresh OldMagneticField from the current face B (the CT update computes
       MagneticField = OldMagneticField - Curl(E)); EvolveLevel does this
       before each hydro solve. */
    g.CopyBaryonFieldToOldBaryonField();
    rc = g.SolveHydroEquations(0, 0, NULL, 0);
    double mab = g.EnzoModulesMaxDivB();
    if (mab > maxafter) maxafter = mab;
  }
  *maxdivb_after = maxafter;

  if (rc != FAIL) {
    g.EnzoModulesGetField(iDen, d);
    g.EnzoModulesGetField(iTE,  e);
    g.EnzoModulesGetField(iV1,  u);
    g.EnzoModulesGetField(iV2,  v);
    g.EnzoModulesGetField(iV3,  w);
    g.EnzoModulesGetField(iB1,  bx);
    g.EnzoModulesGetField(iB2,  by);
    g.EnzoModulesGetField(iB3,  bz);
  }
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
