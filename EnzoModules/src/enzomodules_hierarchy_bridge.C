/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- AMR hierarchy inter-grid operators
/
/  PURPOSE:
/    Exercise the grid:: methods that tie AMR refinement levels together,
/    in isolation, via the generic grid fixture (Grid_EnzoModulesFixture.C):
/
/      1. RESTRICTION  -- grid::ProjectSolutionToParentGrid:
/         volume-average a fine child field down onto its coarse parent
/         (conservation: each parent cell = mean of the RefineBy^rank child
/         cells it covers).
/
/      2. PROLONGATION -- grid::InterpolateFieldValues:
/         interpolate a coarse parent field up to the fine child (the same
/         conservative/monotone interpolator used to fill new subgrids;
/         2nd-order-accurate, so exact for a linear field).
/
/      3. REFLUXING    -- grid::CorrectForRefinedFluxes:
/         conservative flux correction of the coarse cells adjacent to a
/         fine/coarse boundary, given Initial vs. Refined boundary fluxes.
/
/    A coarse PARENT grid covers the unit domain at resolution N; a fine
/    CHILD grid covers the parent active cells [N/4, 3N/4) at 2x resolution
/    (CellWidth_child = CellWidth_parent / 2), with the child edges aligned
/    to parent cell boundaries.  The methods derive the parent<->child index
/    mapping from the edges + CellWidth, so getting the edges right is the
/    crux; RefineBy and the interpolation globals are set explicitly.
/
/    Compiled with Enzo's headers, linked against the Enzo shared library.
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

/* Configure the AMR/hydro globals these inter-grid operators read.  Refinement
 * factor 2, Cartesian, PPM (non-Zeus so density-weighted conservative project /
 * interpolate / reflux paths are taken), SecondOrderA interpolation (exact for
 * linear fields), conservative interpolation on.  Both fixture grids live on
 * processor 0 with SEND_RECEIVE so CommunicationMethodShouldExit proceeds. */
static void em_hier_globals(TopGridData &MetaData)
{
  SetDefaultGlobalValues(MetaData);
  enzomodules_init_timer();

  ComovingCoordinates       = 0;
  NumberOfGhostZones        = 3;
  RefineBy                  = 2;
  Coordinate                = Cartesian;
  HydroMethod               = PPM_DirectEuler;
  UseHydro                  = 1;
  UseMHDCT                  = 0;
  DualEnergyFormalism       = 0;
  RadiativeCooling          = 0;
  FluxCorrection            = 1;
  Gamma                     = 1.4;
  InterpolationMethod       = SecondOrderA;
  ConservativeInterpolation = 1;
  CorrectParentBoundaryFlux = FALSE;
  MyProcessorNumber         = 0;
  NumberOfProcessors        = 1;
  CommunicationDirection    = COMMUNICATION_SEND_RECEIVE;

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainLeftEdge[dim]  = 0.0;
    DomainRightEdge[dim] = 1.0;
  }
}

/* Build the parent + child pair sharing the unit domain.  The parent has `n`
 * active cells over [0,1]; the child covers parent active cells [n/4, 3n/4) at
 * 2x resolution.  Both grids carry Density,TotalEnergy,Velocity1-3 so the
 * conservative (density-weighted) code paths and IdentifyPhysicalQuantities
 * succeed.  Returns parent/child total (incl. ghost) cell counts. */
static const int NHFIELDS = 5;
static void em_build_pair(grid &parent, grid &child, int n,
                          int *parent_size, int *child_size,
                          int *child_active)
{
  const int ng = NumberOfGhostZones;
  int ftypes[NHFIELDS] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };

  /* parent: active [0,1], dims = n + 2*ng */
  int   pdims[3]  = { n + 2 * ng, 1, 1 };
  FLOAT pleft[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT pright[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  *parent_size = parent.EnzoModulesSetupGrid(1, pdims, pleft, pright,
                                             NHFIELDS, ftypes, 0.0);

  /* child: covers parent active cells [n/4, 3n/4) -> nf = (n/2)*2 active cells
     at half the parent cell width. */
  double dxp = 1.0 / (double)n;
  int    nf  = (n / 2) * 2;     /* = n; child active-cell count */
  *child_active = nf;
  double cleft  = (n / 4) * dxp;
  double cright = cleft + nf * (dxp / 2.0);

  int   cdims[3]  = { nf + 2 * ng, 1, 1 };
  FLOAT clo[3]    = { (FLOAT)cleft,  (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT chi[3]    = { (FLOAT)cright, (FLOAT)1.0, (FLOAT)1.0 };
  *child_size = child.EnzoModulesSetupGrid(1, cdims, clo, chi,
                                           NHFIELDS, ftypes, 0.0);
}

extern "C" {

/* ---- RESTRICTION: grid::ProjectSolutionToParentGrid ------------------- */

/* Fill the CHILD density with caller-supplied values (length child_size, incl.
 * ghosts), project it down to the PARENT, and return the PARENT density
 * (length parent_size) plus the index mapping needed to certify conservation:
 *   *out_pstart = parent active-cell index of the first overlapped cell,
 *   *out_n      = number of overlapped parent cells (= child_active / 2),
 *   *out_ng     = NumberOfGhostZones.
 * Each overlapped parent cell should equal the mean of the 2 child cells it
 * covers.  Returns 0 on success. */
int enzomodules_project_to_parent(int n, const double *child_density,
                                  double *parent_density_out, int parent_cap,
                                  int *out_pstart, int *out_n, int *out_ng)
{
  TopGridData MetaData;
  em_hier_globals(MetaData);

  grid parent, child;
  int psize, csize, cactive;
  em_build_pair(parent, child, n, &psize, &csize, &cactive);
  if (psize > parent_cap) return 2;

  int icd = child.EnzoModulesFieldIndex(Density);
  child.EnzoModulesSetField(icd, child_density);

  if (child.ProjectSolutionToParentGrid(parent) == FAIL) return 1;

  int ipd = parent.EnzoModulesFieldIndex(Density);
  parent.EnzoModulesGetField(ipd, parent_density_out);

  *out_pstart = n / 4;             /* parent active index of first overlap */
  *out_n      = cactive / RefineBy;
  *out_ng     = NumberOfGhostZones;
  return 0;
}

/* ---- PROLONGATION: grid::InterpolateFieldValues ----------------------- */

/* Fill the PARENT density with caller-supplied values (length parent_size,
 * incl. ghosts), interpolate it up to the entire CHILD grid, and return the
 * CHILD density (length child_size).  *out_ng / *out_cactive let the caller
 * locate the child active cells and reconstruct child cell centers; the parent
 * left edge of the overlap is parent active index n/4, child dx = (1/n)/2.
 * Returns 0 on success. */
int enzomodules_interpolate_to_child(int n, const double *parent_density,
                                     double *child_density_out, int child_cap,
                                     int *out_cactive, int *out_ng,
                                     double *out_cleft, double *out_cdx)
{
  TopGridData MetaData;
  em_hier_globals(MetaData);

  grid parent, child;
  int psize, csize, cactive;
  em_build_pair(parent, child, n, &psize, &csize, &cactive);
  if (csize > child_cap) return 2;

  int ipd = parent.EnzoModulesFieldIndex(Density);
  parent.EnzoModulesSetField(ipd, parent_density);
  /* the interpolator multiplies field*dens then divides by dens for the
     conservative companion fields, so dens must be > 0; density itself is
     interpolated directly.  The fixture zeroes the companion fields, which is
     fine for a pure density-accuracy check. */

  double dxp = 1.0 / (double)n;
  if (child.InterpolateFieldValues(&parent, NULL, &MetaData) == FAIL) return 1;

  int icd = child.EnzoModulesFieldIndex(Density);
  child.EnzoModulesGetField(icd, child_density_out);

  *out_cactive = cactive;
  *out_ng      = NumberOfGhostZones;
  *out_cleft   = (n / 4) * dxp;     /* child active left edge (problem coords) */
  *out_cdx     = dxp / 2.0;
  return 0;
}

/* ---- REFLUXING: grid::CorrectForRefinedFluxes ------------------------- */

/* Exercise the conservative flux correction at a fine/coarse boundary in 1D.
 *
 * We build a single coarse grid (the "parent") of `n` active cells and a fluxes
 * pair describing a subgrid that spans parent active cells [lo, hi] (global
 * indices on the parent level, measured from DomainLeftEdge).  The InitialFlux
 * is what the parent computed; the RefinedFlux is the (better) flux from the
 * fine level.  CorrectForRefinedFluxes adjusts the coarse cells just OUTSIDE
 * the subgrid's left/right faces:
 *     cell left  of subgrid:  Density += (InitialLeft  - RefinedLeft)
 *     cell right of subgrid:  Density -= (InitialRight - RefinedRight)
 *
 * Inputs: density_in (length grid_size incl. ghosts), the subgrid global index
 * range [lo,hi], and the four flux values (init/refined x left/right) for the
 * Density field.  All other-field fluxes are set equal init==refined (no-op).
 * Returns the corrected density in density_out, plus the active indices of the
 * two corrected coarse cells in *out_left_cell / *out_right_cell.  Returns 0 on
 * success. */
int enzomodules_correct_refined_fluxes(int n, const double *density_in,
                                       int lo, int hi,
                                       double init_left, double refined_left,
                                       double init_right, double refined_right,
                                       double *density_out, int grid_cap,
                                       int *out_left_cell, int *out_right_cell)
{
  TopGridData MetaData;
  em_hier_globals(MetaData);
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = reflecting;
  }

  const int ng = NumberOfGhostZones;
  grid g;
  int   gdims[3]  = { n + 2 * ng, 1, 1 };
  FLOAT gleft[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT gright[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  int   ftypes[NHFIELDS] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };
  int gsize = g.EnzoModulesSetupGrid(1, gdims, gleft, gright,
                                     NHFIELDS, ftypes, 1.0);
  if (gsize > grid_cap) return 2;

  /* density: caller-supplied; companion fields kept positive/quiescent so the
     conservative (mult-by-density / divide-by-density) round trip is benign and
     no positivity "undo" fires. */
  int iDen = g.EnzoModulesFieldIndex(Density);
  int iTE  = g.EnzoModulesFieldIndex(TotalEnergy);
  g.EnzoModulesSetField(iDen, density_in);
  std::vector<double> ones(gsize, 1.0);
  g.EnzoModulesSetField(iTE, ones.data());

  /* Build the three fluxes structs.  In 1D only dim 0 carries data.  The flux
     planes are 1x1x1 (single cell each face).  Global indices are on the parent
     level measured from DomainLeftEdge: the subgrid spans active cells [lo,hi]
     so its left face is at global index lo and its right face at hi+1. */
  fluxes Initial, Refined, Boundary;
  InitializeFluxes(&Initial);
  InitializeFluxes(&Refined);
  InitializeFluxes(&Boundary);

  const int dim = 0;
  for (int j = 0; j < MAX_DIMENSION; j++) {
    Initial.LeftFluxStartGlobalIndex[dim][j]  = 0;
    Initial.LeftFluxEndGlobalIndex[dim][j]    = 0;
    Initial.RightFluxStartGlobalIndex[dim][j] = 0;
    Initial.RightFluxEndGlobalIndex[dim][j]   = 0;
  }
  /* left face at global index lo, right face at global index hi+1 */
  Initial.LeftFluxStartGlobalIndex[dim][dim]  = lo;
  Initial.LeftFluxEndGlobalIndex[dim][dim]    = lo;
  Initial.RightFluxStartGlobalIndex[dim][dim] = hi;
  Initial.RightFluxEndGlobalIndex[dim][dim]   = hi;
  for (int j = 0; j < MAX_DIMENSION; j++) {
    Refined.LeftFluxStartGlobalIndex[dim][j]  = Initial.LeftFluxStartGlobalIndex[dim][j];
    Refined.LeftFluxEndGlobalIndex[dim][j]    = Initial.LeftFluxEndGlobalIndex[dim][j];
    Refined.RightFluxStartGlobalIndex[dim][j] = Initial.RightFluxStartGlobalIndex[dim][j];
    Refined.RightFluxEndGlobalIndex[dim][j]   = Initial.RightFluxEndGlobalIndex[dim][j];
  }

  /* allocate one-cell flux planes for every baryon field on dim 0; default
     init==refined (no-op) for all but Density. */
  const int nfields = NHFIELDS;
  for (int field = 0; field < nfields; field++) {
    Initial.LeftFluxes[field][dim]  = new float[1];
    Initial.RightFluxes[field][dim] = new float[1];
    Refined.LeftFluxes[field][dim]  = new float[1];
    Refined.RightFluxes[field][dim] = new float[1];
    Boundary.LeftFluxes[field][dim]  = new float[1];
    Boundary.RightFluxes[field][dim] = new float[1];
    Initial.LeftFluxes[field][dim][0]  = 0.0;
    Initial.RightFluxes[field][dim][0] = 0.0;
    Refined.LeftFluxes[field][dim][0]  = 0.0;
    Refined.RightFluxes[field][dim][0] = 0.0;
    Boundary.LeftFluxes[field][dim][0]  = 0.0;
    Boundary.RightFluxes[field][dim][0] = 0.0;
  }
  Initial.LeftFluxes[iDen][dim][0]  = (float)init_left;
  Refined.LeftFluxes[iDen][dim][0]  = (float)refined_left;
  Initial.RightFluxes[iDen][dim][0] = (float)init_right;
  Refined.RightFluxes[iDen][dim][0] = (float)refined_right;

  int rc = g.CorrectForRefinedFluxes(&Initial, &Refined, &Boundary,
                                     /*SUBlingGrid*/ FALSE, &MetaData);

  if (rc != FAIL)
    g.EnzoModulesGetField(iDen, density_out);

  /* CorrectForRefinedFluxes deletes RefinedFluxes' planes itself; free the
     rest. */
  for (int field = 0; field < nfields; field++) {
    delete [] Initial.LeftFluxes[field][dim];
    delete [] Initial.RightFluxes[field][dim];
    delete [] Boundary.LeftFluxes[field][dim];
    delete [] Boundary.RightFluxes[field][dim];
  }

  /* corrected coarse cells: active index just left of the subgrid (lo-1) and
     just right of it (hi+1). */
  *out_left_cell  = lo - 1;
  *out_right_cell = hi + 1;
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
