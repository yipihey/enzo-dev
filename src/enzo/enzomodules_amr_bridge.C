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

/* Flag cells for refinement by density slope.  `density` is the full field
 * (length prod(dims), incl. ghosts); out_flag receives the FlaggingField
 * (1 = flagged), *out_count the number of flagged cells.  Returns 0. */
int enzomodules_flag_cells(int rank, int dims[], double dx,
                           const double *density, double slope_threshold,
                           int *out_flag, int *out_count)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  ComovingCoordinates = 0;
  NumberOfGhostZones  = 3;
  for (int m = 0; m < MAX_FLAGGING_METHODS; m++)
    CellFlaggingMethod[m] = INT_UNDEFINED;
  CellFlaggingMethod[0]        = 1;                 /* refine by slope */
  SlopeFlaggingFields[0]       = INT_UNDEFINED;     /* all fields */
  MinimumSlopeForRefinement[0] = slope_threshold;

  grid g;
  FLOAT left[3], right[3];
  em_edges(rank, dims, dx, left, right);
  int ftypes[1] = { Density };
  g.EnzoModulesSetupGrid(rank, dims, left, right, 1, ftypes, 0.0);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density), density);

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
