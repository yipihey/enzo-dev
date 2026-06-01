/***********************************************************************
/
/  GRID CLASS: ENZOMODULES UNIT-TEST FIXTURE SUPPORT
/
/  PURPOSE:
/    Generic primitives for building a minimal, self-contained grid from flat
/    arrays so a single grid:: solver method (ZEUS, radiation transfer,
/    gravity, ...) can be exercised in isolation by EnzoModules, and for
/    copying field data back out.  This is the grid-method analogue of the
/    libyt ConvertToLibyt bridge: it lets external (ctypes) callers drive an
/    Enzo grid method without a full simulation.
/
/    Nothing numerical lives here -- only grid construction and field I/O,
/    reusing the existing PrepareGrid / AllocateGrids / FindField machinery.
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

int FindField(int field, int farray[], int numfields);

int grid::EnzoModulesGridSize()
{
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  return size;
}

/* Configure dimensions/edges/fields, allocate the baryon fields, and mark the
 * grid local to this processor.  field_types[] are FieldType enum values
 * (Density, TotalEnergy, Velocity1, ...).  Returns the total cell count. */
int grid::EnzoModulesSetupGrid(int rank, int dims[], FLOAT left[], FLOAT right[],
                               int num_fields, int field_types[], double dt)
{
  /* Field list must be set before AllocateGrids(). */
  NumberOfBaryonFields = num_fields;
  for (int f = 0; f < num_fields; f++)
    FieldType[f] = field_types[f];

  /* Dimensions, edges, ghost-aware active range, CellWidth (PrepareGrid sets
     dims 0..rank-1; force the unused trailing dimensions to a single zone). */
  PrepareGrid(rank, dims, left, right, 0);
  for (int dim = rank; dim < MAX_DIMENSION; dim++) {
    GridDimension[dim]  = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim]   = 0;
  }

  AllocateGrids();           /* allocates + zeroes BaryonField[0..n-1] */

  dtFixed = dt;
  Time = 0.0;
  OldTime = 0.0;
  ProcessorNumber = MyProcessorNumber;   /* "owned" by this rank */

  return this->EnzoModulesGridSize();
}

void grid::EnzoModulesSetField(int field_index, const double *data)
{
  int size = this->EnzoModulesGridSize();
  for (int i = 0; i < size; i++)
    BaryonField[field_index][i] = (float)data[i];
}

void grid::EnzoModulesGetField(int field_index, double *data)
{
  int size = this->EnzoModulesGridSize();
  for (int i = 0; i < size; i++)
    data[i] = (double)BaryonField[field_index][i];
}

int grid::EnzoModulesFieldIndex(int field_type)
{
  return FindField(field_type, FieldType, NumberOfBaryonFields);
}
