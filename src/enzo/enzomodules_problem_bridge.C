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
int InitializeNew(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData,
                  ExternalBoundary &Exterior, float *Initialdt);
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
#endif
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

void enzomodules_problem_get_field(void *h, int gi, int fi, double *out)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetField(fi, out); }

int enzomodules_problem_num_particles(void *h, int gi)
{ return ((EMProblem *)h)->grids[gi]->ReturnNumberOfParticles(); }

void enzomodules_problem_get_particle_pos(void *h, int gi, int dim, double *out)
{ ((EMProblem *)h)->grids[gi]->EnzoModulesGetParticlePosition(dim, out); }

/* ---- Persistent session: drive the timestep loop step-by-step ---------
 * These let a host language (e.g. a Python re-implementation of EvolveLevel)
 * own the time loop and call each orchestration step on the LIVE hierarchy,
 * rather than handing the whole run to EvolveHierarchy.  The handle is the same
 * EMProblem (use the field accessors above to read the evolving state). */

/* Initialize a problem and build its LevelArray for stepping.  Returns a
 * handle, or NULL on failure. */
void *enzomodules_session_init(const char *paramfile)
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

  char *fname = strdup(paramfile);
  int rc = InitializeNew(fname, p->TopGrid, p->MetaData, p->Exterior, &p->dt);
  if (rc == FAIL) { free(fname); delete p; return NULL; }

  p->MetaData.dtDataDump = 0.0;
  /* The root-grid Poisson FFT defaults to the MPI-only transpose; force the
   * serial path since this library is built without MPI. */
  UnigridTranspose = 0;
  AddLevel(p->LevelArray, &p->TopGrid, 0);

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
 * The radiation sources come from RadiativeTransferInitialize (which read them
 * from the parameter file in session_init).  We deliberately do NOT call
 * RadiativeTransferPrepare here: its StarParticleRadTransfer rebuilds the source
 * list from star particles and would wipe the parameter-file sources.  Instead
 * we size dtPhoton with RadiativeTransferComputeTimestep (which leaves the
 * source list untouched) and let EvolvePhotons sub-cycle and emit. */
int enzomodules_session_evolve_photons(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
#ifdef TRANSFER
  /* RestartPhotons (gated on this flag) would also rebuild from star particles;
   * clear it so EvolvePhotons just transports the existing sources. */
  p->MetaData.FirstTimestepAfterRestart = FALSE;

  /* Build the per-cell SubgridMarker (grid-ownership map the photon transport
   * uses for grid-to-grid handoff).  Normally done inside RebuildHierarchy;
   * the photon walk dereferences it, so it must exist. */
  SetSubgridMarker(p->MetaData, p->LevelArray, level, FALSE);

  /* Size the photon timestep (light-crossing / CFL based) without disturbing
   * the source list. */
  RadiativeTransferComputeTimestep(p->LevelArray, &p->MetaData, 0.0, level);

  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  FLOAT GridTime = (n > 0)
    ? Grids[0]->GridData->ReturnTime() + Grids[0]->GridData->ReturnTimeStep()
    : p->MetaData.Time;
  delete[] Grids;

  Star *AllStars = NULL;
  /* LoopTime=1: sub-cycle PhotonTime up to GridTime in dtPhoton steps so the
   * sources emit (the first sub-cycle just advances PhotonTime past the source
   * creation time; subsequent ones deposit). */
  int rc = EvolvePhotons(&p->MetaData, p->LevelArray, AllStars, GridTime, level, 1);
  return (rc == FAIL) ? 1 : 0;
#else
  return 0;
#endif
}

/* CFL timestep for `level`: min over grids, clamped to StopTime. */
double enzomodules_session_compute_dt(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  double dt = 1e30;
  for (int i = 0; i < n; i++) {
    double dtg = (double)Grids[i]->GridData->ComputeTimeStep();
    if (dtg < dt) dt = dtg;
  }
  delete[] Grids;
  double remaining = (double)p->MetaData.StopTime - (double)p->MetaData.Time;
  if (dt > remaining && remaining > 0) dt = remaining;
  return dt;
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
int enzomodules_session_solve_hydro(void *h, int level)
{
  EMProblem *p = (EMProblem *)h;
  HierarchyEntry **Grids;
  int n = GenerateGridArray(p->LevelArray, level, &Grids);
  int *nsub = p->NumberOfSubgrids[level];
  fluxes ***flux = p->SubgridFluxesEstimate[level];
  int rc = 0;
  for (int i = 0; i < n; i++) {
    int ns = nsub ? nsub[i] : 0;
    fluxes **sf = flux ? flux[i] : NULL;
    if (Grids[i]->GridData->SolveHydroEquations(p->MetaData.CycleNumber, ns,
                                                sf, level) == FAIL) rc = 1;
  }
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
