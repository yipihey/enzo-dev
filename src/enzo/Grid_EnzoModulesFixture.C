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
FLOAT FindCrossSection(int type, float energy);

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

/* ---- Particles -------------------------------------------------------- */

int grid::EnzoModulesSetupParticles(int n, int num_attributes)
{
  NumberOfParticles = n;
  NumberOfParticleAttributes = num_attributes;   /* global */
  this->AllocateNewParticles(n);
  for (int i = 0; i < n; i++) {
    ParticleMass[i]   = 0.0;
    ParticleNumber[i] = i;
    ParticleType[i]   = PARTICLE_TYPE_DARK_MATTER;
    for (int dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim][i] = 0.0;
      ParticleVelocity[dim][i] = 0.0;
    }
    for (int a = 0; a < num_attributes; a++)
      ParticleAttribute[a][i] = 0.0;
  }
  return n;
}

void grid::EnzoModulesSetParticlePosition(int dim, const double *data)
{
  for (int i = 0; i < NumberOfParticles; i++)
    ParticlePosition[dim][i] = (FLOAT)data[i];
}

void grid::EnzoModulesSetParticleVelocity(int dim, const double *data)
{
  for (int i = 0; i < NumberOfParticles; i++)
    ParticleVelocity[dim][i] = (float)data[i];
}

void grid::EnzoModulesSetParticleMass(const double *data)
{
  for (int i = 0; i < NumberOfParticles; i++)
    ParticleMass[i] = (float)data[i];
}

void grid::EnzoModulesGetParticlePosition(int dim, double *data)
{
  for (int i = 0; i < NumberOfParticles; i++)
    data[i] = (double)ParticlePosition[dim][i];
}

/* CIC-deposit the grid's particles onto its own GravitatingMassFieldParticles
 * (the standard particle-mesh operation), and expose the result for tests. */
int grid::EnzoModulesDepositParticles()
{
  GravityBoundaryType = TopGridPeriodic;   /* deposit field at grid resolution */
  this->InitializeGravitatingMassFieldParticles(1);
  this->ClearGravitatingMassFieldParticles();
  if (this->DepositParticlePositions(this, Time,
                                     GRAVITATING_MASS_FIELD_PARTICLES) == FAIL)
    return -1;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];
  return size;
}

void grid::EnzoModulesGetDepositField(double *data)
{
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GravitatingMassFieldParticlesDimension[dim];
  for (int i = 0; i < size; i++)
    data[i] = (double)GravitatingMassFieldParticles[i];
}

double grid::EnzoModulesDepositCellVolume()
{
  double vol = 1.0;
  for (int dim = 0; dim < GridRank; dim++)
    vol *= (double)GravitatingMassFieldParticlesCellSize;
  return vol;
}

/* ---- Radiative transfer (single photon-package ray) ------------------ */

int grid::EnzoModulesRaytrace(double energy, double photons, double dtphoton,
                              double lightspeed, int ipix, int hpix_level,
                              double *photons_final, double *radius_final,
                              double *kph_sum)
{
  int size = this->EnzoModulesGridSize();

  /* Periodic so a single-grid ray wraps and stays in the uniform medium. */
  GravityBoundaryType = TopGridPeriodic;

  /* Every cell of this single grid is owned by this grid. */
  SubgridMarker = new grid *[size];
  for (int i = 0; i < size; i++) SubgridMarker[i] = this;
  HasRadiation = FALSE;
  MaximumkphIfront = 0.0;
  IndexOfMaximumkph = 0;

  /* A sentinel head plus one real package, properly linked (the walk checks
     PreviousPackage->NextPackage == package). */
  PhotonPackages = new PhotonPackageEntry;
  PhotonPackages->NextPackage = NULL;
  PhotonPackages->PreviousPackage = NULL;
  PhotonPackages->CurrentSource = NULL;

  PhotonPackageEntry *P = new PhotonPackageEntry;
  P->PreviousPackage = PhotonPackages;
  P->NextPackage     = NULL;
  PhotonPackages->NextPackage = P;
  P->Photons   = photons;
  P->Type      = 0;                  /* HI-ionizing */
  P->Energy    = energy;
  P->CrossSection = FindCrossSection(0, (float)energy);   /* cm^2 */
  P->EmissionTimeInterval = dtphoton;
  P->EmissionTime = 0.0;
  P->CurrentTime  = 0.0;
  P->Radius       = 0.0;
  P->ColumnDensity = 0.0;
  P->ipix  = ipix;
  P->level = hpix_level;
  P->CurrentSource = NULL;
  for (int dim = 0; dim < 3; dim++)
    P->SourcePosition[dim] = 0.5 * (GridLeftEdge[dim] + GridRightEdge[dim]);
  P->SourcePositionDiff = 0.0;

  grid  *MoveToGrid = NULL;
  grid  *Grids0[1]  = { this };
  int    DeleteMe = FALSE, PauseMe = FALSE, DeltaLevel = 0;
  float  LightCrossingTime = (float)(2.0 * dtphoton);

  PhotonPackageEntry *PP = P;
  this->WalkPhotonPackage(&PP, &MoveToGrid, NULL, this, Grids0, 1,
                          DeleteMe, PauseMe, DeltaLevel,
                          LightCrossingTime, (float)lightspeed,
                          hpix_level, /*MinimumPhotonFlux*/ 0.0);

  *photons_final = (double)P->Photons;
  *radius_final  = (double)P->Radius;

  int kphHINum = FindField(kphHI, FieldType, NumberOfBaryonFields);
  double s = 0.0;
  if (kphHINum >= 0)
    for (int i = 0; i < size; i++) s += (double)BaryonField[kphHINum][i];
  *kph_sum = s;

  delete P;
  delete PhotonPackages;
  delete[] SubgridMarker;
  SubgridMarker = NULL;
  PhotonPackages = NULL;
  return 0;
}
