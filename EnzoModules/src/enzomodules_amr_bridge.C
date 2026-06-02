/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- AMR control (cell flagging + clustering)
/
/  PURPOSE:
/    Exercise the adaptive-mesh-refinement control path in isolation:
/      1. cell flagging  -- grid::SetFlaggingField marks cells needing
/         refinement (here, by density slope);
/      2. clustering      -- the Berger-Rigoutsos signature algorithm
/         (ProtoSubgrid + IdentifyNewSubgridsBySignature) groups flagged cells
/         into minimal rectangular subgrids (the "generate new grids" step --
/         each ProtoSubgrid defines a child grid).
/
/    Compiled with Enzo headers, linked against libenzo.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
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
#include "ProtoSubgrid.h"

int SetDefaultGlobalValues(TopGridData &MetaData);
int IdentifyNewSubgridsBySignature(ProtoSubgrid *SubgridList[],
                                   int &NumberOfSubgrids);

static void em_edges(int rank, int dims[], double dx, FLOAT left[], FLOAT right[])
{
  for (int d = 0; d < 3; d++) {
    left[d]  = 0.0;
    right[d] = (d < rank) ? (FLOAT)((dims[d] - 2 * NumberOfGhostZones) * dx)
                          : (FLOAT)1.0;
  }
}

extern "C" {

/* Flag cells for refinement by a chosen CellFlaggingMethod criterion, via the
 * real grid::SetFlaggingField dispatch.  Supported `method` codes (Enzo's):
 *   1  slope            2  baryon mass / overdensity   3  shocks
 *   9  shear            15 second derivative           (others reachable too,
 *   but need particles/metals/B/chemistry -- set those fields and call the
 *   matching method).  `threshold` sets that criterion's parameter.  d/e/u/v/w
 *   are the full hydro fields (length prod(dims), incl. ghosts); out_flag gets
 *   the FlaggingField (>=1 = flagged), *out_count the flagged-cell count. */
int enzomodules_flag_cells(int rank, int dims[], double dx, int method,
                           double threshold,
                           const double *d, const double *e,
                           const double *u, const double *v, const double *w,
                           int *out_flag, int *out_count)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;
  DualEnergyFormalism = 0;
  Gamma               = 1.4;
  for (int m = 0; m < MAX_FLAGGING_METHODS; m++)
    CellFlaggingMethod[m] = INT_UNDEFINED;
  CellFlaggingMethod[0] = method;

  switch (method) {
    case 1:   /* slope */
      SlopeFlaggingFields[0]       = INT_UNDEFINED;
      MinimumSlopeForRefinement[0] = threshold;
      break;
    case 2:   /* baryon mass / overdensity */
      MinimumMassForRefinement[0]              = threshold;
      MinimumMassForRefinementLevelExponent[0] = 0.0;
      break;
    case 3:   /* shocks */
      MinimumPressureJumpForRefinement = threshold;
      break;
    case 9:   /* shear */
      MinimumShearForRefinement = threshold;
      break;
    case 15:  /* second derivative */
      SecondDerivativeFlaggingFields[0]       = INT_UNDEFINED;
      MinimumSecondDerivativeForRefinement[0] = threshold;
      break;
    default:
      break;
  }

  grid g;
  FLOAT left[3], right[3];
  em_edges(rank, dims, dx, left, right);
  int ftypes[5] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };
  g.EnzoModulesSetupGrid(rank, dims, left, right, 5, ftypes, 0.0);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density),     d);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(TotalEnergy), e);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity1),   u);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity2),   v);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity3),   w);

  g.ClearFlaggingField();
  int count = 0;
  g.SetFlaggingField(count, 0);
  g.EnzoModulesGetFlagging(out_flag);
  *out_count = count;
  return 0;
}

/* Flag cells in a must-refine region (geometric criterion, method 12).
 * region_left/right are the box edges in domain coordinates. */
int enzomodules_flag_region(int rank, int dims[], double dx,
                            const double *region_left,
                            const double *region_right,
                            int *out_flag, int *out_count)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;
  for (int m = 0; m < MAX_FLAGGING_METHODS; m++)
    CellFlaggingMethod[m] = INT_UNDEFINED;
  CellFlaggingMethod[0] = 12;                 /* must-refine region */
  MustRefineRegionMinRefinementLevel = 1;     /* > our level (0) */
  for (int d = 0; d < MAX_DIMENSION; d++) {
    MustRefineRegionLeftEdge[d]  = (d < rank) ? (FLOAT)region_left[d]  : (FLOAT)0.0;
    MustRefineRegionRightEdge[d] = (d < rank) ? (FLOAT)region_right[d] : (FLOAT)1.0;
  }

  grid g;
  FLOAT left[3], right[3];
  em_edges(rank, dims, dx, left, right);
  int ftypes[1] = { Density };
  g.EnzoModulesSetupGrid(rank, dims, left, right, 1, ftypes, 0.0);
  std::vector<double> ones(g.EnzoModulesGridSize(), 1.0);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density), ones.data());

  g.ClearFlaggingField();
  int count = 0;
  g.SetFlaggingField(count, 0);
  g.EnzoModulesGetFlagging(out_flag);
  *out_count = count;
  return 0;
}

/* Flag cells by Jeans length (method 6): refine where the Jeans length is
 * under-resolved (dense/cold gas).  units are physical cgs scalings. */
int enzomodules_flag_jeans(int rank, int dims[], double dx,
                           const double *density, const double *energy,
                           double safety, double density_units,
                           double length_units, double time_units,
                           int *out_flag, int *out_count)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;
  DualEnergyFormalism = 0;
  Gamma               = 5.0 / 3.0;
  Mu                  = 0.6;
  MultiSpecies        = 0;
  RefineByJeansLengthSafetyFactor = safety;
  GlobalDensityUnits = density_units;
  GlobalLengthUnits  = length_units;
  GlobalTimeUnits    = time_units;
  GlobalMassUnits    = density_units * length_units * length_units * length_units;
  for (int m = 0; m < MAX_FLAGGING_METHODS; m++)
    CellFlaggingMethod[m] = INT_UNDEFINED;
  CellFlaggingMethod[0] = 6;                  /* refine by Jeans length */

  grid g;
  FLOAT left[3], right[3];
  em_edges(rank, dims, dx, left, right);
  int ftypes[5] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3 };
  g.EnzoModulesSetupGrid(rank, dims, left, right, 5, ftypes, 0.0);
  std::vector<double> zero(g.EnzoModulesGridSize(), 0.0);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density),     density);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(TotalEnergy), energy);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity1),   zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity2),   zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity3),   zero.data());

  g.ClearFlaggingField();
  int count = 0;
  g.SetFlaggingField(count, 0);
  g.EnzoModulesGetFlagging(out_flag);
  *out_count = count;
  return 0;
}

/* Cluster a flagging field into rectangular subgrids (Berger-Rigoutsos).
 * `flagging` is the input FlaggingField (length prod(dims)); returns the
 * number of subgrids in *out_nsub and, for each, its [lo,hi] active-cell index
 * box per dimension in out_boxes (laid out [sub][2*dim + {0,1}], 6 ints/sub).
 * Returns 0. */
int enzomodules_cluster(int rank, int dims[], double dx, const int *flagging,
                        int *out_nsub, int *out_boxes, int max_sub)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;

  grid g;
  FLOAT left[3], right[3];
  em_edges(rank, dims, dx, left, right);
  int ftypes[1] = { Density };
  g.EnzoModulesSetupGrid(rank, dims, left, right, 1, ftypes, 0.0);
  g.ClearFlaggingField();
  g.EnzoModulesSetFlagging(flagging);

  ProtoSubgrid *SubgridList[MAX_NUMBER_OF_SUBGRIDS];
  int NumberOfSubgrids = 1;
  SubgridList[0] = new ProtoSubgrid;
  SubgridList[0]->SetLevel(1);
  if (SubgridList[0]->CopyFlaggedZonesFromGrid(&g) == FAIL) return 1;
  if (IdentifyNewSubgridsBySignature(SubgridList, NumberOfSubgrids) == FAIL)
    return 2;

  *out_nsub = NumberOfSubgrids;
  for (int s = 0; s < NumberOfSubgrids && s < max_sub; s++) {
    FLOAT *le = SubgridList[s]->ReturnGridLeftEdge();
    FLOAT *re = SubgridList[s]->ReturnGridRightEdge();
    for (int d = 0; d < 3; d++) {
      /* active-cell index range covered by this subgrid in the parent grid */
      int lo = (d < rank) ? (int)(((double)le[d]) / dx + 0.5) : 0;
      int hi = (d < rank) ? (int)(((double)re[d]) / dx + 0.5) : 1;
      out_boxes[s * 6 + 2 * d + 0] = lo;
      out_boxes[s * 6 + 2 * d + 1] = hi;
    }
  }
  for (int s = 0; s < NumberOfSubgrids; s++) delete SubgridList[s];
  return 0;
}

} /* extern "C" */
