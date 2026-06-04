/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- problem setup & initial conditions
/
/  PURPOSE:
/    Drive Enzo's problem initialization (InitializeNew) from a parameter
/    file and expose the resulting grid hierarchy (fields, particles) to
/    external callers.  Because InitializeNew dispatches on ProblemType to
/    every problem generator, this single bridge can build the initial
/    conditions of *any* Enzo problem type from its .enzo parameter file.
/
/    A stateful handle holds the initialized hierarchy; accessors walk it and
/    copy out field / particle data using the grid fixture primitives.
/
/    Compiled with Enzo's headers and linked against the Enzo shared library.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "FastSiblingLocator.h"

class ImplicitProblemABC;
extern "C" void enzomodules_init_timer();
int CommunicationInitialize(Eint32 *argc, char **argv[]);
#ifdef USE_MPI
Eflt64 CommunicationMinValue(Eflt64 Value);   /* global min reduction (CommunicationUtilities.C) */
int CommunicationAllSumValues(Eflt64 *Values, int Number);   /* global sum to all ranks */
#endif
int InitializeNew(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData,
                  ExternalBoundary &Exterior, float *Initialdt);
class ImplicitProblemABC;
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
                       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
                       ImplicitProblemABC *ImplicitSolver,
#endif
                       FLOAT WriteTime, int CheckpointDump);
int Group_ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
                      ExternalBoundary *Exterior, FLOAT *Initialdt,
                      bool ReadParticlesOnly);
int SetDefaultGlobalValues(TopGridData &MetaData);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
                      HierarchyEntry **Grids[]);
int CreateSiblingList(HierarchyEntry **Grids, int NumberOfGrids,
                      SiblingGridList *SiblingList, int StaticLevelZero,
                      TopGridData *MetaData, int level);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          SiblingGridList SiblingList[], int level,
                          TopGridData *MetaData, ExternalBoundary *Exterior,
                          LevelHierarchyEntry *Level);
int RebuildHierarchy(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                     int level);
int PrepareDensityField(LevelHierarchyEntry *LevelArray[], int level,
                        TopGridData *MetaData, FLOAT When,
                        SiblingGridList **SiblingGridListStorage);
int UpdateParticlePositions(grid *Grid);
int CreateFluxes(HierarchyEntry *Grids[], fluxes **SubgridFluxesEstimate[],
                 int NumberOfGrids, int NumberOfSubgrids[]);
int FinalizeFluxes(HierarchyEntry *Grids[], fluxes **SubgridFluxesEstimate[],
                   int NumberOfGrids, int NumberOfSubgrids[]);
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
                         int NumberOfSubgrids[],
                         fluxes **SubgridFluxesEstimate[],
                         LevelHierarchyEntry *SUBlingList[],
                         TopGridData *MetaData);
int CreateSUBlingList(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                      int level, SiblingGridList SiblingList[],
                      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids, LevelHierarchyEntry **SUBlingList);
class Star;
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                  Star *&AllStars, FLOAT GridTime, int level, int LoopTime);
#ifdef TRANSFER
int RadiativeTransferInitialize(char *ParameterFile, HierarchyEntry &TopGrid,
                                TopGridData &MetaData, ExternalBoundary &Exterior,
                                ImplicitProblemABC *&ImplicitSolver,
                                LevelHierarchyEntry *LevelArray[]);
int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
                                     TopGridData *MetaData, float dtLevelAbove,
                                     int level);
int SetSubgridMarker(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
                     int level, int UpdateReplicatedGridsOnly);
int StarParticleRadTransfer(LevelHierarchyEntry *LevelArray[], int level,
                            Star *AllStars);
#endif
/* Star formation / feedback / activation lifecycle (not TRANSFER-specific). */
int StarParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
                           int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
                           int ThisLevel, Star *&AllStars,
                           int TotalStarParticleCountPrevious[]);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
                         int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
                         int level, Star *&AllStars,
                         int TotalStarParticleCountPrevious[], int &OutputNow);
int ComputeDednerWaveSpeeds(TopGridData *MetaData,
                            LevelHierarchyEntry *LevelArray[], int level,
                            double dt0);
int ActiveParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
                             int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
                             int ThisLevel);
int ActiveParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
                           int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
                           int level, int NumberOfNewActiveParticles[]);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
                         TopGridData *MetaData);
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float *norm, float *pTopGridTimeStep);
int ComputeStochasticForcing(TopGridData *MetaData, HierarchyEntry *Grids[],
                             int NumberOfGrids);
int ComputeDomainBoundaryMassFlux(HierarchyEntry *Grids[], int level,
                                  int NumberOfGrids, TopGridData *MetaData);
int CallProblemSpecificRoutines(TopGridData *MetaData, HierarchyEntry *ThisGrid,
                                int GridNum, float *norm, float TopGridTimeStep,
                                int level, int LevelCycleCount[]);
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
#ifdef TRANSFER
                    ImplicitProblemABC *ImplicitSolver,
#endif
                    LevelHierarchyEntry *LevelArray[], float Initialdt);

struct EMProblem {
  HierarchyEntry TopGrid;
  TopGridData    MetaData;
  ExternalBoundary Exterior;
  float dt;
  std::vector<grid *> grids;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  SiblingGridList *SiblingGridListStorage[MAX_DEPTH_OF_HIERARCHY];
  /* Per-level transient flux state, threaded CreateFluxes -> SolveHydroEquations
   * -> UpdateFromFinerGrids -> FinalizeFluxes (the AMR conservation machinery). */
  int *NumberOfSubgrids[MAX_DEPTH_OF_HIERARCHY];
  fluxes ***SubgridFluxesEstimate[MAX_DEPTH_OF_HIERARCHY];
};

/* Depth-first walk of the AMR hierarchy collecting every grid on this rank. */
static void collect_grids(HierarchyEntry *he, std::vector<grid *> &out)
{
  while (he != NULL) {
    if (he->GridData != NULL)
      out.push_back(he->GridData);
    if (he->NextGridNextLevel != NULL)
      collect_grids(he->NextGridNextLevel, out);
    he = he->NextGridThisLevel;
  }
}

extern "C" {

/* Initialize a problem from its parameter file.  Returns an opaque handle, or
 * NULL on failure.  Output files (OutputLog, ...) are written to the current
 * working directory, so callers should run from a scratch directory. */
void *enzomodules_init_problem(const char *paramfile)
{
  static bool comm_done = false;
  if (!comm_done) {
    int argc = 1;
    static char arg0[] = "enzomodules";
    static char *argv_storage[2] = { arg0, NULL };
    char **argv = argv_storage;
    CommunicationInitialize(&argc, &argv);
    comm_done = true;
  }

  EMProblem *p = new EMProblem();
  p->TopGrid.NextGridThisLevel = NULL;
  p->TopGrid.NextGridNextLevel = NULL;
  p->TopGrid.ParentGrid        = NULL;
  p->TopGrid.GridData          = NULL;

  char *fname = strdup(paramfile);
  int rc = InitializeNew(fname, p->TopGrid, p->MetaData, p->Exterior, &p->dt);
  free(fname);
  if (rc == FAIL) { delete p; return NULL; }

  collect_grids(&p->TopGrid, p->grids);
  return (void *)p;
}

/* Initialize a problem, then run Enzo's full time integration
 * (EvolveHierarchy: per-grid timesteps, boundary conditions, the hydro/MHD/
 * gravity solvers, AMR sub-cycling, flux correction and projection) until
 * StopTime/StopCycle.  stop_time>0 / stop_cycle>0 override the parameter file.
 * Data dumps are suppressed.  Returns the evolved-hierarchy handle (read with
 * the same accessors), or NULL on failure. */
void *enzomodules_evolve_problem(const char *paramfile, double stop_time,
                                 int stop_cycle)
{
  static bool comm_done = false;
  if (!comm_done) {
    int argc = 1;
    static char arg0[] = "enzomodules";
    static char *argv_storage[2] = { arg0, NULL };
    char **argv = argv_storage;
    CommunicationInitialize(&argc, &argv);
    comm_done = true;
  }

  EMProblem *p = new EMProblem();
  p->TopGrid.NextGridThisLevel = NULL;
  p->TopGrid.NextGridNextLevel = NULL;
  p->TopGrid.ParentGrid        = NULL;
  p->TopGrid.GridData          = NULL;

  char *fname = strdup(paramfile);
  int rc = InitializeNew(fname, p->TopGrid, p->MetaData, p->Exterior, &p->dt);
  free(fname);
  if (rc == FAIL) { delete p; return NULL; }

  enzomodules_init_timer();              /* EvolveHierarchy uses TIMER_START */
  if (stop_time  > 0.0) p->MetaData.StopTime  = stop_time;
  if (stop_cycle > 0)   p->MetaData.StopCycle = stop_cycle;
  p->MetaData.dtDataDump        = 0.0;   /* 0 = never (suppress output) */
  p->MetaData.CycleSkipDataDump = 0;

  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  for (int l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) LevelArray[l] = NULL;
  AddLevel(LevelArray, &p->TopGrid, 0);

  if (EvolveHierarchy(p->TopGrid, p->MetaData, &p->Exterior,
#ifdef TRANSFER
                      NULL,
#endif
                      LevelArray, p->dt) == FAIL) {
    delete p;
    return NULL;
  }

  collect_grids(&p->TopGrid, p->grids);
  return (void *)p;
}

int enzomodules_problem_problemtype(void *)       { return ProblemType; }
int enzomodules_problem_num_grids(void *h)         { return (int)((EMProblem *)h)->grids.size(); }

/* The refinement level of grid `gi` (0 = root), or -1 if not in the hierarchy.
 * EMProblem->grids and LevelArray reference the same grid objects (both via
 * collect_grids over the hierarchy), so a per-level pointer-membership scan
 * recovers the level — the enabling primitive for a :julia hydro slot to iterate
 * the grids on a given level under AMR. */
int enzomodules_problem_grid_level(void *h, int gi)
{
  EMProblem *p = (EMProblem *)h;
  grid *target = p->grids[gi];
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    HierarchyEntry **Grids;
    int n = GenerateGridArray(p->LevelArray, level, &Grids);
    int found = -1;
    for (int i = 0; i < n; i++)
      if (Grids[i]->GridData == target) { found = level; break; }
    delete[] Grids;
    if (found >= 0) return found;
    if (n == 0) break;                 /* no grids at this level ⇒ none deeper */
  }
  return -1;
}

int enzomodules_problem_grid_rank(void *h, int gi)
{ return ((EMProblem *)h)->grids[gi]->GetGridRank(); }

void enzomodules_problem_grid_dims(void *h, int gi, int *dims)
{
  grid *g = ((EMProblem *)h)->grids[gi];
  for (int d = 0; d < 3; d++) dims[d] = g->GetGridDimension(d);
}

int enzomodules_problem_num_fields(void *h, int gi)
{ return ((EMProblem *)h)->grids[gi]->ReturnNumberOfBaryonFields(); }

void enzomodules_problem_field_types(void *h, int gi, int *types)
{ ((EMProblem *)h)->grids[gi]->ReturnFieldType(types); }

int enzomodules_problem_grid_size(void *h, int gi)
{ return ((EMProblem *)h)->grids[gi]->EnzoModulesGridSize(); }

/* MPI locality: this rank's id and the rank count (0 / 1 in the serial flavor),
 * and the home processor of grid `gi` (always 0 in serial).  These let the Julia
 * side iterate only the grids resident on this rank, exactly as Enzo's own
 * SolveHydroEquations does (it skips grids whose ProcessorNumber != mine). */
int enzomodules_session_my_rank(void *h)   { return MyProcessorNumber; }
int enzomodules_session_num_ranks(void *h) { return NumberOfProcessors; }

int enzomodules_problem_grid_processor(void *h, int gi)
{
  EMProblem *p = (EMProblem *)h;
  if (gi < 0 || gi >= (int)p->grids.size()) return -1;
  return p->grids[gi]->ReturnProcessorNumber();
}

void enzomodules_problem_get_field(void *h, int gi, int fi, double *out)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetField(fi, out); }

/* Write a field (flat, incl. ghost zones) back into the LIVE grid's BaryonField.
 * Enables a host-language physics method to mutate the Enzo state in place (the
 * ':julia' slot swap): read with get_field, compute, write back here. */
void enzomodules_problem_set_field(void *h, int gi, int fi, const double *in)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesSetField(fi, in); }

/* Self-gravity coupling: write/read the cell-centered AccelerationField[dim], so a
 * :julia gravity slot can compute g = -grad(phi) and feed Enzo's SolveHydroEquations. */
void enzomodules_problem_set_acceleration(void *h, int gi, int dim, const double *in)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesSetAcceleration(dim, in); }

void enzomodules_problem_get_acceleration(void *h, int gi, int dim, double *out)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetAcceleration(dim, out); }

/* ---- ADR-0003 part B: conservative :julia hydro under AMR ---------------
 * Bridge for writing EnzoNG's recorded face fluxes into Enzo's flux registers,
 * so a :julia hydro slot can feed UpdateFromFinerGrids/CorrectForRefinedFluxes:
 *   - each grid's BoundaryFluxes  = the RefinedFluxes a finer grid carried, and
 *   - the parent's SubgridFluxesEstimate[level][i][sub] = the coarse InitialFluxes
 *     under subgrid `sub`.
 * Together these are exactly the two flux sets SolveHydroEquations fills; with
 * them, Enzo's own machinery restores conservation across coarse-fine boundaries. */

/* The i-th grid on `level` in GenerateGridArray order (the SAME order
 * create_fluxes used to build SubgridFluxesEstimate[level]). */
static grid *em_grid_on_level(EMProblem *p, int level, int i)
{
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  grid *g = (i >= 0 && i < n) ? Grids[i]->GridData : NULL;
  delete[] Grids;
  return g;
}

/* Flat grid-list index (the gi the field/mesh accessors use) of the i-th grid on
 * `level`, so a :julia AMR slot can build an EnzoGridMesh for the same grid whose
 * subgrid fluxes it fills by (level, i). -1 if out of range. */
int enzomodules_problem_grid_index_on_level(void *h, int level, int i)
{
  EMProblem *p = (EMProblem *)h;
  grid *g = em_grid_on_level(p, level, i);
  if (g == NULL) return -1;
  for (size_t k = 0; k < p->grids.size(); k++)
    if (p->grids[k] == g) return (int)k;
  return -1;
}

/* Grid geometry / BoundaryFluxes accessors (grid-method wrappers, by flat gi). */
void enzomodules_problem_grid_global_start(void *h, int gi, long_int *gstart)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGlobalStart(gstart); }

void enzomodules_problem_grid_edge(void *h, int gi, double *left, double *right)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGridEdge(left, right); }

int enzomodules_problem_boundary_flux_size(void *h, int gi, int dim)
{ return ((EMProblem *)h)->grids[gi]->EnzoModulesBoundaryFluxSize(dim); }

void enzomodules_problem_boundary_flux_extent(void *h, int gi, int dim, int side,
                                              long_int *start, long_int *end)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesBoundaryFluxExtent(dim, side, start, end); }

void enzomodules_problem_set_boundary_flux(void *h, int gi, int field, int dim,
                                           int side, const double *plane)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesSetBoundaryFlux(field, dim, side, plane); }

void enzomodules_problem_get_boundary_flux(void *h, int gi, int field, int dim,
                                           int side, double *plane)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetBoundaryFlux(field, dim, side, plane); }

/* Number of subgrid flux entries for the i-th grid on `level` (proper subgrids
 * + 1 for the grid's own external boundary, the last entry). Requires
 * create_fluxes(level) first; -1 if absent. */
int enzomodules_problem_num_subgrids(void *h, int level, int i)
{
  EMProblem *p = (EMProblem *)h;
  if (p->NumberOfSubgrids[level] == NULL) return -1;
  return p->NumberOfSubgrids[level][i];
}

/* Coarse-index global extents of subgrid flux (level,i,sub) for `dim`/`side`
 * (the footprint ReturnFluxDims set from the child's BoundaryFluxes / refinement).
 * Length-3 outputs. */
void enzomodules_problem_subgrid_flux_extent(void *h, int level, int i, int sub,
                                             int dim, int side,
                                             long_int *start, long_int *end)
{
  EMProblem *p = (EMProblem *)h;
  fluxes *F = p->SubgridFluxesEstimate[level][i][sub];
  long_int *s = (side == 0) ? F->LeftFluxStartGlobalIndex[dim]
                            : F->RightFluxStartGlobalIndex[dim];
  long_int *e = (side == 0) ? F->LeftFluxEndGlobalIndex[dim]
                            : F->RightFluxEndGlobalIndex[dim];
  for (int d = 0; d < 3; d++) { start[d] = s[d]; end[d] = e[d]; }
}

/* Size (plane cells) of subgrid flux (level,i,sub) for `dim` — the coarse-index
 * footprint extent (1 in 1D). */
int enzomodules_problem_subgrid_flux_size(void *h, int level, int i, int sub, int dim)
{
  EMProblem *p = (EMProblem *)h;
  grid *g = em_grid_on_level(p, level, i);
  int rank = g ? g->GetGridRank() : 1;
  fluxes *F = p->SubgridFluxesEstimate[level][i][sub];
  int size = 1;
  for (int d = 0; d < rank; d++)
    size *= F->LeftFluxEndGlobalIndex[dim][d] - F->LeftFluxStartGlobalIndex[dim][d] + 1;
  return size;
}

/* SET (overwrite) a coarse InitialFlux plane into SubgridFluxesEstimate
 * [level][i][sub]->{Left,Right}Fluxes[field][dim]. Allocates the array if NULL
 * (CreateFluxes/ReturnFluxDims leaves the flux pointers NULL; SolveHydroEquations
 * normally allocates them — we do the same). Overwrite (not add): the coarse
 * InitialFlux is one coarse step's flux, recreated each step by create_fluxes. */
void enzomodules_problem_set_subgrid_flux(void *h, int level, int i, int sub,
                                          int field, int dim, int side,
                                          const double *plane)
{
  EMProblem *p = (EMProblem *)h;
  int size = enzomodules_problem_subgrid_flux_size(h, level, i, sub, dim);
  fluxes *F = p->SubgridFluxesEstimate[level][i][sub];
  float **arr = (side == 0) ? F->LeftFluxes[field] : F->RightFluxes[field];
  if (arr[dim] == NULL) arr[dim] = new float[size];
  for (int n = 0; n < size; n++) arr[dim][n] = (float)plane[n];
}

void enzomodules_problem_get_subgrid_flux(void *h, int level, int i, int sub,
                                          int field, int dim, int side, double *plane)
{
  EMProblem *p = (EMProblem *)h;
  int size = enzomodules_problem_subgrid_flux_size(h, level, i, sub, dim);
  fluxes *F = p->SubgridFluxesEstimate[level][i][sub];
  float **arr = (side == 0) ? F->LeftFluxes[field] : F->RightFluxes[field];
  for (int n = 0; n < size; n++) plane[n] = (arr[dim] == NULL) ? 0.0 : (double)arr[dim][n];
}

int enzomodules_problem_num_particles(void *h, int gi)
{ return ((EMProblem *)h)->grids[gi]->ReturnNumberOfParticles(); }

void enzomodules_problem_get_particle_pos(void *h, int gi, int dim, double *out)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetParticlePosition(dim, out); }

/* ---- Persistent session: drive the timestep loop step-by-step ---------
 * These let a host language (e.g. a Python re-implementation of EvolveLevel)
 * own the time loop and call each orchestration step on the LIVE hierarchy,
 * rather than handing the whole run to EvolveHierarchy.  The handle is the same
 * EMProblem (use the field accessors above to read the evolving state). */

/* Internal accessors so sibling bridge translation units (e.g. the inline
 * halo-finder bridge) can reach the live hierarchy on a session handle without
 * duplicating the EMProblem layout.  extern "C" to dodge name mangling; the C++
 * return types just need to be declared identically in the caller. */
extern "C" LevelHierarchyEntry **EnzoModulesProblemLevelArray(void *h)
{ return ((EMProblem *)h)->LevelArray; }

extern "C" TopGridData *EnzoModulesProblemMetaData(void *h)
{ return &((EMProblem *)h)->MetaData; }


/* One-time MPI/comm + timer init, shared by session_init and from_output. */
static void em_ensure_comm()
{
  static bool comm_done = false;
  if (!comm_done) {
    int argc = 1;
    static char arg0[] = "enzomodules";
    static char *argv_storage[2] = { arg0, NULL };
    char **argv = argv_storage;
    CommunicationInitialize(&argc, &argv);
    comm_done = true;
  }
  enzomodules_init_timer();
}

/* Allocate an EMProblem with its hierarchy/level state zeroed. */
static EMProblem *em_new_problem()
{
  EMProblem *p = new EMProblem();
  p->TopGrid.NextGridThisLevel = NULL;
  p->TopGrid.NextGridNextLevel = NULL;
  p->TopGrid.ParentGrid        = NULL;
  p->TopGrid.GridData          = NULL;
  for (int l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) {
    p->LevelArray[l] = NULL;
    p->SiblingGridListStorage[l] = NULL;
    p->NumberOfSubgrids[l] = NULL;
    p->SubgridFluxesEstimate[l] = NULL;
  }
  return p;
}

/* Initialize a problem and build its LevelArray for stepping.  Returns a
 * handle, or NULL on failure. */
void *enzomodules_session_init(const char *paramfile)
{
  em_ensure_comm();

  EMProblem *p = em_new_problem();

  /* Hardening: several radiative-transfer globals are reset only by
   * RadiativeTransferReadParameters (the RT-on path), not by
   * SetDefaultGlobalValues.  Because this library keeps one set of Enzo globals
   * for the whole process, they leak from a prior RT session into a later
   * non-RT one -- e.g. RadiationPressure left at 1 makes the PPM solver expect a
   * gravity AccelerationField that a pure-hydro problem never built, and it
   * segfaults.  Reset them up front so each session starts clean; an RT problem
   * re-reads its own values in RadiativeTransferInitialize below. */
#ifdef TRANSFER
  RadiationPressure = FALSE;
#endif

  char *fname = strdup(paramfile);
  int rc;
  /* InitializeNew throws EnzoFatalException (ENZO_FAIL) on many setup errors --
   * e.g. a missing/unparseable cooling-rate data file.  Catch it so a failed
   * problem init returns NULL gracefully instead of aborting the host process. */
  try {
    rc = InitializeNew(fname, p->TopGrid, p->MetaData, p->Exterior, &p->dt);
  } catch (EnzoFatalException &e) {
    rc = FAIL;
  }
  if (rc == FAIL) { free(fname); delete p; return NULL; }

  p->MetaData.dtDataDump = 0.0;
#ifndef USE_MPI
  /* The root-grid Poisson FFT defaults to the MPI-only transpose; force the
   * serial path in the serial flavor of this library.  In the MPI flavor we
   * keep the parameter-file/transpose default so the distributed FFT runs. */
  UnigridTranspose = 0;
#endif
  AddLevel(p->LevelArray, &p->TopGrid, 0);

  /* Initial-grid distribution across ranks is done by InitializeNew's partition
   * loop (CommunicationPartitionGrid, a no-op when NumberOfProcessors==1), and
   * regrid load-balancing by RebuildHierarchy — both already invoked on the
   * paths below, so no extra distribution call is needed here. */

#ifdef TRANSFER
  /* Radiative transfer needs its own init (done in enzo.C, not InitializeNew):
   * it allocates the kph/PhotoGamma fields, sets dtPhoton, and builds the
   * radiation source list.  Skipped cleanly when RadiativeTransfer is off. */
  if (RadiativeTransfer) {
    ImplicitProblemABC *ImplicitSolver = NULL;
    RadiativeTransferInitialize(fname, p->TopGrid, p->MetaData, p->Exterior,
                                ImplicitSolver, p->LevelArray);
  }
#endif
  free(fname);

  collect_grids(&p->TopGrid, p->grids);
  return (void *)p;
}

/* Write the full simulation state to disk (Group_WriteAllData) -- the same HDF5
 * data dump / checkpoint enzo.C writes, with grid data, the hierarchy, the
 * external boundary and a parameter file.  `basename` + zero-padded `filenumber`
 * name the dump (e.g. basename="snap", filenumber=0 -> "snap0000").  Set
 * checkpoint != 0 for a checkpoint dump (flags CheckpointRestart in the output).
 * Returns 0 on success, 1 on failure.  The resulting dump can be reloaded with
 * enzomodules_session_from_output. */
int enzomodules_session_write_output(void *h, int filenumber, int checkpoint)
{
  EMProblem *p = (EMProblem *)h;
  /* Group_WriteAllData builds its output name as <name-buffer> + <id(filenumber)>,
   * but only initializes the name buffer when the basename matches one of the
   * recognized dump-name patterns -- otherwise it is uninitialized stack memory
   * (garbage prefix).  Pass the session's own DataDumpName so the DataDumpName
   * branch fires, and NULL the dump directories so the dump lands directly in
   * the working directory (no subdir) with the deterministic name
   * "<DataDumpName><id>".  Restore the fields afterwards. */
  const char *base = (p->MetaData.DataDumpName && p->MetaData.DataDumpName[0])
                   ? p->MetaData.DataDumpName : "data";
  char *name = strdup(base);
  char *saved_ddir = p->MetaData.DataDumpDir;
  char *saved_gdir = p->MetaData.GlobalDir;
  char *saved_ldir = p->MetaData.LocalDir;
  p->MetaData.DataDumpDir = NULL;
  p->MetaData.GlobalDir   = NULL;
  p->MetaData.LocalDir    = NULL;
  int rc;
  try {
    rc = Group_WriteAllData(name, filenumber, &p->TopGrid, p->MetaData,
                            &p->Exterior,
#ifdef TRANSFER
                            NULL,   /* ImplicitSolver: FLD only, unused for Moray RT */
#endif
                            -1, checkpoint ? TRUE : FALSE);
  } catch (EnzoFatalException &e) {
    rc = FAIL;
  }
  p->MetaData.DataDumpDir = saved_ddir;
  p->MetaData.GlobalDir   = saved_gdir;
  p->MetaData.LocalDir    = saved_ldir;
  free(name);
  return (rc == FAIL) ? 1 : 0;
}

/* Reload a simulation from a dump written by enzomodules_session_write_output
 * (or by Enzo itself): Group_ReadAllData reads the parameter file, external
 * boundary, hierarchy and grid data, then we build the LevelArray (all levels)
 * exactly as session_init does.  `name` is the dump's hierarchy/parameter file
 * (e.g. "snap0000").  Returns a session handle, or NULL on failure. */
void *enzomodules_session_from_output(const char *name)
{
  em_ensure_comm();

  EMProblem *p = em_new_problem();

  /* Reset globals first (Group_ReadAllData's ReadParameterFile only sets what
   * the dump lists), then clear the RT global that leaks into non-RT runs. */
  SetDefaultGlobalValues(p->MetaData);
#ifdef TRANSFER
  RadiationPressure = FALSE;
#endif
  UnigridTranspose = 0;

  char *fname = strdup(name);
  FLOAT initdt = 0.0;
  int rc;
  try {
    rc = Group_ReadAllData(fname, &p->TopGrid, p->MetaData, &p->Exterior,
                           &initdt, false);
  } catch (EnzoFatalException &e) {
    rc = FAIL;
  }
  if (rc == FAIL) { free(fname); delete p; return NULL; }
  p->dt = (float)initdt;

  p->MetaData.dtDataDump = 0.0;
  AddLevel(p->LevelArray, &p->TopGrid, 0);

#ifdef TRANSFER
  if (RadiativeTransfer) {
    ImplicitProblemABC *ImplicitSolver = NULL;
    RadiativeTransferInitialize(fname, p->TopGrid, p->MetaData, p->Exterior,
                                ImplicitSolver, p->LevelArray);
  }
#endif
  free(fname);

  collect_grids(&p->TopGrid, p->grids);
  return (void *)p;
}

double enzomodules_session_time(void *h)      { return (double)((EMProblem *)h)->MetaData.Time; }
double enzomodules_session_stop_time(void *h) { return (double)((EMProblem *)h)->MetaData.StopTime; }
int    enzomodules_session_cycle(void *h)     { return ((EMProblem *)h)->MetaData.CycleNumber; }

/* Refresh the cached grid list (call after rebuild). */
static void em_recollect(EMProblem *p)
{
  p->grids.clear();
  collect_grids(&p->TopGrid, p->grids);
}

/* Apply boundary conditions (sibling/parent copies + external boundary) to all
 * grids on `level`. */
int enzomodules_session_set_boundary(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  /* libenzo is built with -DFAST_SIB, so SetBoundaryConditions takes the
   * precomputed sibling list (same path EvolveLevel uses).  Cache it in the
   * per-level storage so the gravity chain (PrepareDensityField) can reuse it. */
  delete[] p->SiblingGridListStorage[level];
  SiblingGridList *SiblingList = new SiblingGridList[n];
  CreateSiblingList(Grids, n, SiblingList, 0, &p->MetaData, level);
  p->SiblingGridListStorage[level] = SiblingList;
  int rc = SetBoundaryConditions(Grids, n, SiblingList, level, &p->MetaData,
                                 &p->Exterior, p->LevelArray[level]);
  delete[] Grids;
  return (rc == FAIL) ? 1 : 0;
}

/* Self-gravity chain for `level` (the same sequence EvolveLevel runs):
 * deposit mass + solve the Poisson equation (PrepareDensityField), then per
 * grid compute accelerations and copy the potential to the baryon field.
 * Requires set_boundary(level) first (to build the sibling list).  No-op
 * unless SelfGravity is on. */
int enzomodules_session_gravity(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  if (p->SiblingGridListStorage[level] == NULL) {
    /* gravity needs a sibling list; build one if set_boundary wasn't called */
    enzomodules_session_set_boundary(h, level);
  }
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  if (PrepareDensityField(p->LevelArray, level, &p->MetaData, 0.5,
                          p->SiblingGridListStorage) == FAIL) rc = 1;
  for (int i = 0; i < n; i++) {
    grid *g = Grids[i]->GridData;
    if (level > 0) g->SolveForPotential(level);
    g->ComputeAccelerations(level);
    g->CopyPotentialToBaryonField();
    g->ComputeAccelerationFieldExternal();
  }
  delete[] Grids;
  return rc;
}

/* Update particle positions on all grids of `level` (drift by dtFixed). */
void enzomodules_session_update_particles(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  for (int i = 0; i < n; i++)
    UpdateParticlePositions(Grids[i]->GridData);
  delete[] Grids;
}

/* Radiative transfer: trace photon packages on `level`, sub-cycling the photon
 * time up to the grid time (GridTime = grid Time + dtFixed) so the sources
 * actually emit and deposit photo-ionization/heating rates.  No-op unless
 * RadiativeTransfer is on (and the library was built with -DTRANSFER).
 *
 * Two source modes:
 *  - use_star_sources == 0 (default): use the radiation sources that
 *    RadiativeTransferInitialize read from the parameter file (e.g. PhotonTest).
 *    We deliberately do NOT call RadiativeTransferPrepare, whose
 *    StarParticleRadTransfer would rebuild the list from star particles and wipe
 *    those sources.
 *  - use_star_sources != 0: gather the grid's star particles into the AllStars
 *    list (StarParticleInitialize) and convert them into radiation sources
 *    (StarParticleRadTransfer), so stars radiate (e.g. ProblemType 252).
 * In both cases we size dtPhoton with RadiativeTransferComputeTimestep and let
 * EvolvePhotons sub-cycle and emit. */
int enzomodules_session_evolve_photons_ex(void *h, int level,
                                          int use_star_sources)
{
  EMProblem *p = (EMProblem *)h;
#ifdef TRANSFER
  /* Hardening: nothing to do if radiative transfer is off or the level holds no
   * grids -- return success rather than walking null state. */
  if (!RadiativeTransfer) return 0;
  if (p->LevelArray[level] == NULL) return 0;

  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  if (n == 0) { delete[] Grids; return 0; }
  Star *AllStars = NULL;

  if (use_star_sources) {
    /* Build the master star list from the grid particles.  FindAllStarParticles
     * (which turns star particles into Star objects) is gated on
     * FirstTimestepAfterRestart, so leave that flag set until after this call. */
    int *TotalPrev = new int[n > 0 ? n : 1];
    for (int i = 0; i < n; i++) TotalPrev[i] = 0;
    StarParticleInitialize(Grids, &p->MetaData, n, p->LevelArray, level,
                           AllStars, TotalPrev);
    delete[] TotalPrev;
  }

  /* RestartPhotons (gated on this flag) would rebuild the field from scratch;
   * clear it so EvolvePhotons just transports the current sources. */
  p->MetaData.FirstTimestepAfterRestart = FALSE;

  /* Build the per-cell SubgridMarker (grid-ownership map the photon transport
   * uses for grid-to-grid handoff).  Normally done inside RebuildHierarchy;
   * the photon walk dereferences it, so it must exist. */
  SetSubgridMarker(p->MetaData, p->LevelArray, level, FALSE);

  /* Size the photon timestep (light-crossing / CFL based). */
  RadiativeTransferComputeTimestep(p->LevelArray, &p->MetaData, 0.0, level);

  if (use_star_sources)
    /* Convert the star particles into the radiation-source list. */
    StarParticleRadTransfer(p->LevelArray, level, AllStars);

  FLOAT GridTime = (n > 0)
    ? Grids[0]->GridData->ReturnTime() + Grids[0]->GridData->ReturnTimeStep()
    : p->MetaData.Time;
  delete[] Grids;

  /* LoopTime=1: sub-cycle PhotonTime up to GridTime in dtPhoton steps so the
   * sources emit (the first sub-cycle just advances PhotonTime past the source
   * creation time; subsequent ones deposit). */
  int rc = EvolvePhotons(&p->MetaData, p->LevelArray, AllStars, GridTime, level, 1);
  return (rc == FAIL) ? 1 : 0;
#else
  return 0;
#endif
}

/* Backward-compatible entry: parameter-file sources. */
int enzomodules_session_evolve_photons(void *h, int level)
{
  return enzomodules_session_evolve_photons_ex(h, level, 0);
}

/* CFL timestep for `level`: min over grids, clamped to StopTime. */
double enzomodules_session_compute_dt(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  double dt = 1e30;
  for (int i = 0; i < n; i++) {
    /* Skip grids owned by another rank: their data is not resident here, and
     * a non-local grid contributes its dt via that rank's reduction below. */
    if (Grids[i]->GridData->ReturnProcessorNumber() != MyProcessorNumber)
      continue;
    double dtg = (double)Grids[i]->GridData->ComputeTimeStep();
    if (dtg < dt) dt = dtg;
  }
  delete[] Grids;
#ifdef USE_MPI
  /* The timestep is a global quantity: reduce the per-rank minimum across all
   * ranks (no-op when NumberOfProcessors==1). */
  dt = (double)CommunicationMinValue((Eflt64)dt);
#endif
  double remaining = (double)p->MetaData.StopTime - (double)p->MetaData.Time;
  if (dt > remaining && remaining > 0) dt = remaining;
  return dt;
}

/* Global composite integral of BaryonField[field] over the ROOT grid (level 0),
 * weighted by cell volume — the conserved total the reflux gate checks (mass for
 * field=Density).  Each rank sums its LOCAL level-0 tiles (the root grid is split
 * one tile per rank by CommunicationPartitionGrid), then the partial sums are
 * reduced to ALL ranks, so every rank returns the same global value.  This is the
 * multi-rank conservation primitive (ADR-0005 #4); in serial it is just the local
 * sum (the All-sum is a no-op when NumberOfProcessors==1). */
double enzomodules_session_global_field_integral(void *h, int field)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, 0, &Grids);
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    if (Grids[i]->GridData->ReturnProcessorNumber() != MyProcessorNumber)
      continue;   /* non-local tile: counted by its owning rank's partial sum */
    total += Grids[i]->GridData->EnzoModulesActiveFieldIntegral(field);
  }
  delete[] Grids;
#ifdef USE_MPI
  Eflt64 v = (Eflt64)total;
  CommunicationAllSumValues(&v, 1);   /* reduce partial sums to all ranks */
  total = (double)v;
#endif
  return total;
}

/* Set the timestep on all grids of `level`. */
void enzomodules_session_set_dt(void *h, int level, double dt)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  for (int i = 0; i < n; i++) Grids[i]->GridData->SetTimeStep((float)dt);
  delete[] Grids;
}

/* Solve the hydro/MHD equations on all grids of `level` (one step).  If
 * create_fluxes(level) was called first, the boundary fluxes are accumulated
 * into the per-level flux storage (needed for conservative AMR flux
 * correction); otherwise no subgrid fluxes are recorded (fine for unigrid). */
/* Advance the hydro/MHD equations one step on `level`.  Dispatches on
 * HydroMethod, mirroring EvolveLevel:
 *  - PPM (0), Zeus (2), constrained-transport MHD (MHD_Li, 6): the single-call
 *    Grid::SolveHydroEquations.
 *  - Runge-Kutta HD (HD_RK, 3) and MHD (MHD_RK, 4): the 2nd-order two-step
 *    integration (1st step, refresh boundaries, 2nd step), with Dedner wave
 *    speeds computed up front for MHD_RK.  (The optional RK2 gravity re-deposit
 *    EvolveLevel does for self-gravitating RK runs is omitted; call gravity()
 *    around solve_hydro for that case.) */
/* Apply the comoving expansion (Hubble drag) source terms on a level's grids.
 * In EvolveLevel this runs once per grid after advance_time when
 * ComovingCoordinates is on (EvolveLevel.C); the step-by-step session loop must
 * call it too for cosmology runs, since it is NOT part of SolveHydroEquations. */
int enzomodules_session_comoving_expansion(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->ComovingExpansionTerms() == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* CT-MHD (UseMHDCT): allocate + zero the per-grid AvgElectricField accumulator.
 * EvolveLevel.C:377 does this for level>0 grids at level entry; the subgrid's
 * SolveHydroEquations sums its electric field into it, and the parent's
 * UpdateFromFinerGrids/MHD_ProjectFace reads it (null without this -> segfault).
 * ClearAvgElectricField is itself a no-op when UseMHDCT is off. */
int enzomodules_session_clear_avg_electric_field(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->ClearAvgElectricField() == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* CT-MHD: after UpdateFromFinerGrids, recompute the cell-centered B from the
 * face-corrected magnetic field using the averaged electric field from finer
 * grids (the CT analog of flux refluxing) -- EvolveLevel.C:899-901, guarded by
 * UseMHDCT && MHD_ProjectE. NextLevel = LevelArray[level+1]. */
int enzomodules_session_mhd_update_magnetic_field(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->MHD_UpdateMagneticField(level, p->LevelArray[level+1],
                                                    FALSE) == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

int enzomodules_session_solve_hydro(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int *nsub = p->NumberOfSubgrids[level];
  fluxes ***flux = p->SubgridFluxesEstimate[level];
  int rc = 0;

  if (HydroMethod != HD_RK && HydroMethod != MHD_RK) {
    for (int i = 0; i < n; i++) {
      int ns = nsub ? nsub[i] : 0;
      fluxes **sf = flux ? flux[i] : NULL;
      if (Grids[i]->GridData->SolveHydroEquations(p->MetaData.CycleNumber, ns,
                                                  sf, level) == FAIL) rc = 1;
    }
    delete[] Grids;
    return rc;
  }

  /* --- Runge-Kutta two-step path (HD_RK / MHD_RK) --- */
  /* Set the global NColor from the grid's fields (default INT_UNDEFINED); the
   * RK solvers size their primitive arrays as NEQ+NSpecies+NColor, so this must
   * run first or those stack arrays get a garbage size.  The RK driver does the
   * same on Grids[0]. */
  if (n > 0) Grids[0]->GridData->SetNumberOfColours();

  double dt = (n > 0) ? (double)Grids[0]->GridData->ReturnTimeStep() : 0.0;
  if (HydroMethod == MHD_RK && level == 0)
    ComputeDednerWaveSpeeds(&p->MetaData, p->LevelArray, level, dt);

  for (int i = 0; i < n; i++) {
    int ns = nsub ? nsub[i] : 0;
    fluxes **sf = flux ? flux[i] : NULL;
    if (HydroMethod == HD_RK) {
      if (Grids[i]->GridData->RungeKutta2_1stStep(sf, ns, level,
                                                  &p->Exterior) == FAIL) rc = 1;
    } else {
      if (Grids[i]->GridData->MHDRK2_1stStep(sf, ns, level,
                                             &p->Exterior) == FAIL) rc = 1;
    }
  }

  /* Refresh boundaries between the RK sub-steps (uses the sibling list built by
   * set_boundary; rebuild it if absent). */
  SiblingGridList *SiblingList = p->SiblingGridListStorage[level];
  if (SiblingList == NULL) {
    SiblingList = new SiblingGridList[n];
    CreateSiblingList(Grids, n, SiblingList, 0, &p->MetaData, level);
    p->SiblingGridListStorage[level] = SiblingList;
  }
  SetBoundaryConditions(Grids, n, SiblingList, level, &p->MetaData,
                        &p->Exterior, p->LevelArray[level]);

  for (int i = 0; i < n; i++) {
    int ns = nsub ? nsub[i] : 0;
    fluxes **sf = flux ? flux[i] : NULL;
    if (HydroMethod == HD_RK) {
      if (Grids[i]->GridData->RungeKutta2_2ndStep(sf, ns, level,
                                                  &p->Exterior) == FAIL) rc = 1;
    } else {
      if (Grids[i]->GridData->MHDRK2_2ndStep(sf, ns, level,
                                             &p->Exterior) == FAIL) rc = 1;
    }
  }
  delete[] Grids;
  return rc;
}

/* Solve the radiative cooling + non-equilibrium species rate equations on
 * `level` (grid::MultiSpeciesHandler) -- the separate cooling/chemistry sub-step
 * EvolveLevel runs after the hydro solve.  Dispatches to Grackle, the coupled
 * rate-and-cool solver, or SolveRateEquations + SolveRadiativeCooling depending
 * on the run's settings.  A clean no-op (returns 0) when both MultiSpecies and
 * RadiativeCooling are off.  Call after solve_hydro(level). */
int enzomodules_session_solve_cooling(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->MultiSpeciesHandler() == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* Run the star-particle lifecycle on `level`, the same sequence EvolveLevel
 * drives: StarParticleInitialize (gather existing stars into the AllStars list +
 * set feedback flags), then per-grid StarParticleHandler (star formation + the
 * feedback deposit, STARMAKE_METHOD), then StarParticleFinalize (ActivateNewStar
 * -- flip newly-formed/aged particles to active radiating stars -- and sync back
 * to the grids).  A clean no-op (returns 0) when both StarParticleCreation and
 * StarParticleFeedback are off.  Call after solve_hydro(level)/solve_cooling so
 * formation sees the updated gas; a later evolve_photons(level, stars=True) then
 * radiates from the now-active stars. */
int enzomodules_session_star_particles(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  if (!StarParticleCreation && !StarParticleFeedback) return 0;

  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  if (n == 0) { delete[] Grids; return 0; }

  Star *AllStars = NULL;
  int *TotalPrev = new int[n];
  for (int i = 0; i < n; i++) TotalPrev[i] = 0;

  /* The grid's own timestep drives formation/feedback; for a single-level
   * session the level dt and the root timestep coincide. */
  float dt = Grids[0]->GridData->ReturnTimeStep();

  int rc = 0;
  /* Gather existing star particles into Star objects + the AllStars list and
   * set their feedback flags (FindAllStarParticles, gated on
   * FirstTimestepAfterRestart, builds the Star objects from the particles). */
  if (StarParticleInitialize(Grids, &p->MetaData, n, p->LevelArray, level,
                             AllStars, TotalPrev) == FAIL) rc = 1;

  /* Star formation + feedback per grid. */
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->StarParticleHandler(Grids[i]->NextGridNextLevel,
                                                level, dt, dt) == FAIL) rc = 1;

  /* Activate newly-formed / aged stars and sync the AllStars changes back to
   * the grid particles (this is where a Pop III particle flips to an active,
   * radiating star). */
  int OutputNow = FALSE;
  if (StarParticleFinalize(Grids, &p->MetaData, n, p->LevelArray, level,
                           AllStars, TotalPrev, OutputNow) == FAIL) rc = 1;

  delete[] TotalPrev;
  delete[] Grids;
  /* Formation/feedback may have created particles and redistributed them. */
  em_recollect(p);
  return rc;
}

/* Active-particle lifecycle on `level` (the modern sink / SmartStar / accretion
 * framework): ActiveParticleInitialize -> per-grid ActiveParticleHandler ->
 * ActiveParticleFinalize.  No-op unless active particles are enabled. */
int enzomodules_session_active_particles(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  if (n == 0) { delete[] Grids; return 0; }

  int *NumNew = new int[n];
  for (int i = 0; i < n; i++) NumNew[i] = 0;
  float dt = Grids[0]->GridData->ReturnTimeStep();

  int rc = 0;
  if (ActiveParticleInitialize(Grids, &p->MetaData, n, p->LevelArray,
                               level) == FAIL) rc = 1;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->ActiveParticleHandler(Grids[i]->NextGridNextLevel,
                                                  level, dt, NumNew[i]) == FAIL)
      rc = 1;
  if (ActiveParticleFinalize(Grids, &p->MetaData, n, p->LevelArray, level,
                             NumNew) == FAIL) rc = 1;
  delete[] NumNew;
  delete[] Grids;
  em_recollect(p);
  return rc;
}

/* Recompute the homogeneous radiation background (RadiationFieldUpdate) on
 * `level` -- the UV background / Haardt-Madau-style field.  No-op unless a
 * RadiationFieldType is active. */
int enzomodules_session_update_radiation_field(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  return (RadiationFieldUpdate(p->LevelArray, level, &p->MetaData) == FAIL)
         ? 1 : 0;
}

/* Turbulence driving: normalize the random forcing (top grid) and, when
 * DrivenFlowProfile is set, compute the stochastic force field via FFT.  Run
 * before solve_hydro so the forcing enters the momentum update. */
int enzomodules_session_random_forcing(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  float norm = 0.0, topdt = (n > 0) ? Grids[0]->GridData->ReturnTimeStep() : 0.0;
  int rc = 0;
  if (ComputeRandomForcingNormalization(p->LevelArray, 0, &p->MetaData,
                                        &norm, &topdt) == FAIL) rc = 1;
  if (DrivenFlowProfile)
    if (ComputeStochasticForcing(&p->MetaData, Grids, n) == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* Per-grid physics steps EvolveLevel runs after the hydro solve: thermal
 * conduction (ConductHeat) and shock finding (ShocksHandler).  No-ops unless the
 * corresponding physics is enabled. */
int enzomodules_session_conduct_heat(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->ConductHeat() == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

int enzomodules_session_find_shocks(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (Grids[i]->GridData->ShocksHandler() == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* Cosmology accessors: compute the expansion factor at the current time.
 * Writes the scale factor a (normalized to 1 at the initial redshift) and the
 * redshift z.  Returns 1 if comoving coordinates are off (a=1, z=0). */
int enzomodules_session_cosmology(void *h, double *a_out, double *z_out)
{
  EMProblem *p = (EMProblem *)h;
  if (!ComovingCoordinates) { if (a_out) *a_out = 1.0;
                              if (z_out) *z_out = 0.0; return 1; }
  FLOAT a = 1.0, dadt = 0.0;
  CosmologyComputeExpansionFactor(p->MetaData.Time, &a, &dadt);
  if (a_out) *a_out = (double)a;
  /* a is in units of (1+InitialRedshift)^-1 ... i.e. a=1 at InitialRedshift, so
   * 1+z = (1+InitialRedshift)/a. */
  if (z_out) *z_out = (double)((1.0 + InitialRedshift) / a - 1.0);
  return 0;
}

/* Track the mass flux through the domain boundary on `level`
 * (ComputeDomainBoundaryMassFlux) -- a conservation diagnostic, accumulated into
 * MetaData's boundary-mass-flux container. */
int enzomodules_session_domain_boundary_mass_flux(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = (ComputeDomainBoundaryMassFlux(Grids, level, n, &p->MetaData) == FAIL)
           ? 1 : 0;
  delete[] Grids;
  return rc;
}

/* Run the problem-specific per-step routines (CallProblemSpecificRoutines) --
 * the ProblemType-dispatched hook EvolveLevel calls each grid step (e.g. the
 * SphericalInfall central reset).  No-op for problem types without one. */
int enzomodules_session_problem_specific_routines(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  float norm = 0.0;
  float topdt = (n > 0) ? Grids[0]->GridData->ReturnTimeStep() : 0.0;
  int rc = 0;
  for (int i = 0; i < n; i++)
    if (CallProblemSpecificRoutines(&p->MetaData, Grids[i], i, &norm, topdt,
                                    level, LevelCycleCount) == FAIL) rc = 1;
  delete[] Grids;
  return rc;
}

/* Allocate the per-level boundary-flux storage (one entry per subgrid) that
 * SolveHydroEquations accumulates into and UpdateFromFinerGrids reads.  Call
 * before solve_hydro(level) when you need conservative AMR flux correction. */
int enzomodules_session_create_fluxes(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  delete[] p->NumberOfSubgrids[level];
  delete[] p->SubgridFluxesEstimate[level];
  p->NumberOfSubgrids[level] = new int[n > 0 ? n : 1];
  p->SubgridFluxesEstimate[level] = new fluxes **[n > 0 ? n : 1];
  int rc = CreateFluxes(Grids, p->SubgridFluxesEstimate[level], n,
                        p->NumberOfSubgrids[level]);
  delete[] Grids;
  return (rc == FAIL) ? 1 : 0;
}

/* Conservative fine->coarse coupling for `level`: project the finer-level
 * solution into these grids and correct their boundary zones for the flux
 * difference (the recursive heart of EvolveLevel).  Requires create_fluxes +
 * solve_hydro on this level and a fully evolved level+1. */
int enzomodules_session_update_from_finer(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  if (p->SiblingGridListStorage[level] == NULL)
    enzomodules_session_set_boundary(h, level);
  LevelHierarchyEntry **SUBlingList = new LevelHierarchyEntry *[n > 0 ? n : 1];
  CreateSUBlingList(&p->MetaData, p->LevelArray, level,
                    p->SiblingGridListStorage[level], &SUBlingList);
  int rc = UpdateFromFinerGrids(level, Grids, n, p->NumberOfSubgrids[level],
                                p->SubgridFluxesEstimate[level], SUBlingList,
                                &p->MetaData);
  DeleteSUBlingList(n, SUBlingList);
  delete[] Grids;
  return (rc == FAIL) ? 1 : 0;
}

/* Release the per-level flux storage allocated by create_fluxes. */
int enzomodules_session_finalize_fluxes(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  if (p->SubgridFluxesEstimate[level] == NULL) return 0;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int rc = FinalizeFluxes(Grids, p->SubgridFluxesEstimate[level], n,
                          p->NumberOfSubgrids[level]);
  delete[] Grids;
  delete[] p->NumberOfSubgrids[level];
  delete[] p->SubgridFluxesEstimate[level];
  p->NumberOfSubgrids[level] = NULL;
  p->SubgridFluxesEstimate[level] = NULL;
  return (rc == FAIL) ? 1 : 0;
}

/* Save the current baryon fields as the "old" fields (time-centering for the
 * gravity source term and the RK schemes). */
void enzomodules_session_copy_baryon_to_old(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  for (int i = 0; i < n; i++)
    Grids[i]->GridData->CopyBaryonFieldToOldBaryonField();
  delete[] Grids;
}

/* Zero the boundary-flux accumulators on all grids of `level` (done once per
 * level step before the subgrid sub-cycles accumulate into them). */
void enzomodules_session_clear_boundary_fluxes(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  for (int i = 0; i < n; i++)
    Grids[i]->GridData->ClearBoundaryFluxes();
  delete[] Grids;
}

/* Number of grids on `level` (so a Python EvolveLevel can iterate/recurse). */
int enzomodules_session_num_grids_on_level(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  delete[] Grids;
  return n;
}

/* Advance each grid's time by its dtFixed and bump the cycle counter. */
void enzomodules_session_advance_time(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  for (int i = 0; i < n; i++) Grids[i]->GridData->SetTimeNextTimestep();
  if (n > 0) p->MetaData.Time = Grids[0]->GridData->ReturnTime();
  p->MetaData.CycleNumber++;
  delete[] Grids;
}

/* Regrid: flag + cluster + (de)refine, then refresh the grid cache. */
int enzomodules_session_rebuild(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  int rc = RebuildHierarchy(&p->MetaData, p->LevelArray, level);
  em_recollect(p);
  return (rc == FAIL) ? 1 : 0;
}

void enzomodules_free_problem(void *h)
{
  EMProblem *p = (EMProblem *)h;
  for (int l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) {
    delete[] p->SiblingGridListStorage[l];
    delete[] p->NumberOfSubgrids[l];
    delete[] p->SubgridFluxesEstimate[l];
  }
  for (size_t i = 0; i < p->grids.size(); i++)
    delete p->grids[i];               /* frees BaryonField / particle arrays */
  delete p;                           /* (subgrid HierarchyEntry nodes leak;
                                         acceptable for a short-lived tool) */
}

} /* extern "C" */
