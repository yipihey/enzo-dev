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
#include <math.h>
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

  /* Hydro-method grid members are normally copied from TopGridData during
     InitializeNew; a bare fixture grid must set them explicitly. */
  PPMFlatteningParameter = 0;
  PPMDiffusionParameter  = 0;
  PPMSteepeningParameter = 0;

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

/* Write the cell-centered AccelerationField[dim] (the gravity source that
 * SolveHydroEquations reads). Allocates the array if it does not exist, since a
 * :julia gravity slot replaces ComputeAccelerations (which normally allocates). */
void grid::EnzoModulesSetAcceleration(int dim, const double *data)
{
  int size = this->EnzoModulesGridSize();
  if (AccelerationField[dim] == NULL)
    AccelerationField[dim] = new float[size];
  for (int i = 0; i < size; i++)
    AccelerationField[dim][i] = (float)data[i];
}

void grid::EnzoModulesGetAcceleration(int dim, double *data)
{
  int size = this->EnzoModulesGridSize();
  for (int i = 0; i < size; i++)
    data[i] = (AccelerationField[dim] == NULL) ? 0.0 : (double)AccelerationField[dim][i];
}

/* ---- ADR-0003 part B: BoundaryFluxes for conservative :julia AMR -------- */

/* This grid's physical edges per dim (length-3 outputs) — so a :julia AMR slot
 * can build a per-grid mesh with the correct (level-dependent) cell width. */
void grid::EnzoModulesGridEdge(double *left, double *right)
{
  for (int dim = 0; dim < GridRank; dim++) {
    left[dim]  = (double)GridLeftEdge[dim];
    right[dim] = (double)GridRightEdge[dim];
  }
  for (int dim = GridRank; dim < MAX_DIMENSION; dim++) { left[dim] = 0.0; right[dim] = 0.0; }
}

/* Global zone index of this grid's first ACTIVE cell per dim (the value
 * PrepareBoundaryFluxes uses for LeftFluxStartGlobalIndex[dim][dim], and the
 * GridGlobalStart-equivalent xEulerSweep offsets the subgrid flux planes by).
 * Lets a :julia hydro map its recorded face fluxes to the global flux indices. */
void grid::EnzoModulesGlobalStart(long_int gstart[])
{
  for (int dim = 0; dim < GridRank; dim++)
    gstart[dim] = nlongint((GridLeftEdge[dim] - DomainLeftEdge[dim]) / CellWidth[dim][0]);
  for (int dim = GridRank; dim < MAX_DIMENSION; dim++)
    gstart[dim] = 0;
}

/* Number of plane cells in this grid's `dim` boundary-flux face (the product of
 * the active extents in the orthogonal dims; 1 in 1D). Ensures the structure
 * exists + is sized/zeroed (ClearBoundaryFluxes). */
int grid::EnzoModulesBoundaryFluxSize(int dim)
{
  if (BoundaryFluxes == NULL || BoundaryFluxes->LeftFluxes[0][dim] == NULL)
    this->ClearBoundaryFluxes();
  int size = 1;
  for (int i = 0; i < GridRank; i++)
    size *= BoundaryFluxes->LeftFluxEndGlobalIndex[dim][i] -
            BoundaryFluxes->LeftFluxStartGlobalIndex[dim][i] + 1;
  return size;
}

/* The global-index extents (start[], end[]) of this grid's `dim`/`side` boundary
 * flux plane (side: 0=Left, 1=Right). Length-3 outputs. */
void grid::EnzoModulesBoundaryFluxExtent(int dim, int side, long_int start[], long_int end[])
{
  if (BoundaryFluxes == NULL)
    this->PrepareBoundaryFluxes();
  long_int *s = (side == 0) ? BoundaryFluxes->LeftFluxStartGlobalIndex[dim]
                            : BoundaryFluxes->RightFluxStartGlobalIndex[dim];
  long_int *e = (side == 0) ? BoundaryFluxes->LeftFluxEndGlobalIndex[dim]
                            : BoundaryFluxes->RightFluxEndGlobalIndex[dim];
  for (int i = 0; i < 3; i++) { start[i] = s[i]; end[i] = e[i]; }
}

/* ADD a boundary-flux plane into this grid's BoundaryFluxes->{Left,Right}Fluxes
 * [field][dim] (accumulate, because a finer grid's outer flux is summed over its
 * temporal subcycles before UpdateFromFinerGrids projects it; ClearBoundaryFluxes
 * zeros once at level entry). Plane length must be EnzoModulesBoundaryFluxSize. */
void grid::EnzoModulesSetBoundaryFlux(int field, int dim, int side, const double *plane)
{
  if (BoundaryFluxes == NULL || BoundaryFluxes->LeftFluxes[field][dim] == NULL)
    this->ClearBoundaryFluxes();
  int size = this->EnzoModulesBoundaryFluxSize(dim);
  float *f = (side == 0) ? BoundaryFluxes->LeftFluxes[field][dim]
                         : BoundaryFluxes->RightFluxes[field][dim];
  for (int i = 0; i < size; i++)
    f[i] += (float)plane[i];
}

void grid::EnzoModulesGetBoundaryFlux(int field, int dim, int side, double *plane)
{
  int size = this->EnzoModulesBoundaryFluxSize(dim);
  float *f = (side == 0) ? BoundaryFluxes->LeftFluxes[field][dim]
                         : BoundaryFluxes->RightFluxes[field][dim];
  for (int i = 0; i < size; i++)
    plane[i] = (f == NULL) ? 0.0 : (double)f[i];
}

int grid::EnzoModulesFieldIndex(int field_type)
{
  return FindField(field_type, FieldType, NumberOfBaryonFields);
}

void grid::EnzoModulesSetFlagging(const int *data)
{
  int size = this->EnzoModulesGridSize();
  for (int i = 0; i < size; i++) FlaggingField[i] = data[i];
}

void grid::EnzoModulesGetFlagging(int *data)
{
  int size = this->EnzoModulesGridSize();
  for (int i = 0; i < size; i++) data[i] = FlaggingField[i];
}

/* ---- MHD constrained transport (face-centered B) --------------------- */

/* Copy a flat row-major array into the face-centered MagneticField[dim].
 * MagneticField[dim] is staggered: it has one extra zone along axis `dim`
 * (MagneticDims[dim] = GridDimension with +1 on axis dim), set up by
 * MHD_SetupDims during EnzoModulesSetupGrid (UseMHDCT must be 1).  The caller
 * supplies a buffer of length MagneticSize[dim]. */
void grid::EnzoModulesSetMagneticField(int dim, const double *data)
{
  for (int i = 0; i < MagneticSize[dim]; i++)
    MagneticField[dim][i] = (float)data[i];
}

int grid::EnzoModulesMagneticSize(int dim)
{
  return MagneticSize[dim];
}

/* Discrete CT divergence of the face-centered magnetic field over the active
 * cells (the defining invariant of constrained transport).  Mirrors the
 * divergence stencil in Grid_MHD_Diagnose.C using the staggered MagneticDims
 * indexing (indexb1/2/3).  Returns the maximum |divB| over active cells. */
double grid::EnzoModulesMaxDivB()
{
  double dx = (double)CellWidth[0][0];
  double dy = (GridRank > 1) ? (double)CellWidth[1][0] : 1.0;
  double dz = (GridRank > 2) ? (double)CellWidth[2][0] : 1.0;

  double maxdivb = 0.0;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
        double divergence =
          (MagneticField[0][indexb1(i+1,j,k)] - MagneticField[0][indexb1(i,j,k)]) / dx +
          ((GridRank < 2) ? 0.0 :
           (MagneticField[1][indexb2(i,j+1,k)] - MagneticField[1][indexb2(i,j,k)]) / dy) +
          ((GridRank < 3) ? 0.0 :
           (MagneticField[2][indexb3(i,j,k+1)] - MagneticField[2][indexb3(i,j,k)]) / dz);
        if (fabs(divergence) > maxdivb) maxdivb = fabs(divergence);
      }
  return maxdivb;
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

/* Multi-grid photon transport: inject a ray into this grid (gA) and let it
 * cross into the sibling gB, mirroring Enzo's inter-grid handoff
 * (SubgridMarker -> FindPhotonNewGrid -> move package to the new grid's list).
 * Verifies AMR photon transport: a ray traversing two tiled grids attenuates
 * exactly as one grid of the combined length. */
int grid::EnzoModulesRaytraceTwoGrid(grid *gB, double energy, double photons,
                                     double dtphoton, double lightspeed,
                                     int ipix, int hpix_level,
                                     double *photons_final, double *radius_final,
                                     int *grids_visited)
{
  int sizeA = this->EnzoModulesGridSize();
  int sizeB = gB->EnzoModulesGridSize();
  FLOAT zero_offset[MAX_DIMENSION] = { 0.0, 0.0, 0.0 };

  this->GravityBoundaryType = SubGridIsolated;
  gB->GravityBoundaryType   = SubGridIsolated;

  /* Every cell -> its own grid, then cells overlapping the sibling -> sibling
     (this is how a ray at the shared face is told to move grids). */
  this->SubgridMarker = new grid *[sizeA];
  for (int i = 0; i < sizeA; i++) this->SubgridMarker[i] = this;
  gB->SubgridMarker = new grid *[sizeB];
  for (int i = 0; i < sizeB; i++) gB->SubgridMarker[i] = gB;
  this->SetSubgridMarkerFromSibling(gB, zero_offset);
  gB->SetSubgridMarkerFromSibling(this, zero_offset);

  this->HasRadiation = FALSE; gB->HasRadiation = FALSE;
  this->MaximumkphIfront = 0.0; this->IndexOfMaximumkph = 0;
  gB->MaximumkphIfront = 0.0;   gB->IndexOfMaximumkph = 0;

  this->PhotonPackages = new PhotonPackageEntry;
  this->PhotonPackages->NextPackage = NULL;
  this->PhotonPackages->PreviousPackage = NULL;
  gB->PhotonPackages = new PhotonPackageEntry;
  gB->PhotonPackages->NextPackage = NULL;
  gB->PhotonPackages->PreviousPackage = NULL;

  PhotonPackageEntry *P = new PhotonPackageEntry;
  P->PreviousPackage = this->PhotonPackages;
  P->NextPackage = NULL;
  this->PhotonPackages->NextPackage = P;
  P->Photons = photons; P->Type = 0; P->Energy = energy;
  P->CrossSection = FindCrossSection(0, (float)energy);
  P->EmissionTimeInterval = dtphoton;
  P->EmissionTime = 0.0; P->CurrentTime = 0.0; P->Radius = 0.0;
  P->ColumnDensity = 0.0; P->ipix = ipix; P->level = hpix_level;
  P->CurrentSource = NULL; P->SourcePositionDiff = 0.0;
  P->SourcePosition[0] = GridLeftEdge[0] + 0.1 * (GridRightEdge[0] - GridLeftEdge[0]);
  P->SourcePosition[1] = 0.5 * (GridLeftEdge[1] + GridRightEdge[1]);
  P->SourcePosition[2] = 0.5 * (GridLeftEdge[2] + GridRightEdge[2]);

  grid *current = this;
  PhotonPackageEntry *PP = P;
  int visited = 1;
  float LightCrossingTime = (float)(4.0 * dtphoton);

  for (int iter = 0; iter < 20; iter++) {
    grid *MoveToGrid = NULL;
    int DeleteMe = FALSE, PauseMe = FALSE, DeltaLevel = 0;
    grid *Grids0[2] = { this, gB };
    current->WalkPhotonPackage(&PP, &MoveToGrid, NULL, current, Grids0, 2,
                               DeleteMe, PauseMe, DeltaLevel,
                               LightCrossingTime, (float)lightspeed,
                               hpix_level, 0.0);
    if (DeleteMe || P->Photons <= 0) break;
    if (MoveToGrid != NULL && MoveToGrid != current) {
      /* relink P out of `current` and into MoveToGrid's list */
      if (P->PreviousPackage) P->PreviousPackage->NextPackage = P->NextPackage;
      if (P->NextPackage) P->NextPackage->PreviousPackage = P->PreviousPackage;
      P->NextPackage = MoveToGrid->PhotonPackages->NextPackage;
      P->PreviousPackage = MoveToGrid->PhotonPackages;
      if (MoveToGrid->PhotonPackages->NextPackage)
        MoveToGrid->PhotonPackages->NextPackage->PreviousPackage = P;
      MoveToGrid->PhotonPackages->NextPackage = P;
      current = MoveToGrid;
      PP = P;
      visited++;
    } else {
      break;
    }
  }

  *photons_final = (double)P->Photons;
  *radius_final  = (double)P->Radius;
  *grids_visited = visited;

  delete P;
  delete this->PhotonPackages; this->PhotonPackages = NULL;
  delete gB->PhotonPackages;   gB->PhotonPackages = NULL;
  delete[] this->SubgridMarker; this->SubgridMarker = NULL;
  delete[] gB->SubgridMarker;   gB->SubgridMarker = NULL;
  return 0;
}
