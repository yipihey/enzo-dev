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

int CommunicationInitialize(Eint32 *argc, char **argv[]);
int InitializeNew(char *filename, HierarchyEntry &TopGrid, TopGridData &MetaData,
                  ExternalBoundary &Exterior, float *Initialdt);

struct EMProblem {
  HierarchyEntry TopGrid;
  TopGridData    MetaData;
  ExternalBoundary Exterior;
  float dt;
  std::vector<grid *> grids;
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

void enzomodules_free_problem(void *h)
{
  EMProblem *p = (EMProblem *)h;
  for (size_t i = 0; i < p->grids.size(); i++)
    delete p->grids[i];               /* frees BaryonField / particle arrays */
  delete p;                           /* (subgrid HierarchyEntry nodes leak;
                                         acceptable for a short-lived tool) */
}

} /* extern "C" */
