/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- INLINE HALO FINDER (FOF) + SUBFIND
/
/  PURPOSE:
/    Expose Enzo's inline friends-of-friends halo finder and the SUBFIND
/    subhalo finder as a certified, in-process session step that a host
/    language (e.g. the Python `Session`) can drive on the LIVE particle
/    hierarchy -- the same FOF/SUBFIND pipeline EvolveHierarchy runs inline,
/    but with its catalogue captured in memory and handed back through getters
/    instead of only written to disk.
/
/    This mirrors the FOF() driver in FOF.C: set units, copy the grid particles
/    into the finder's data structure (FOF_Initialize), build the friends-of-
/    friends groups, optionally run SUBFIND, then read each group's properties
/    with the finder's own get_particles/get_properties (so we reuse the
/    certified property computation rather than re-deriving it).  No behaviour is
/    reimplemented here.
/
/    Single-process only: the MPI slab-exchange / link-across / stitch path in
/    FOF.C is skipped (NumberOfProcessors == 1 takes the same branch there), so
/    this is the serial halo finder.  The caller MUST match the precision this
/    library was built with (see enzomodules_*_precision_bytes()).
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "FOF_allvars.h"
#include "FOF_proto.h"

/* From FOF.C / FOF_Initialize.C / FOF_Finalize.C (compiled into libenzo). */
void FOF_Initialize(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                    FOFData &AllVars, bool SmoothedDarkMatter);
void FOF_Finalize(FOFData &D, LevelHierarchyEntry *LevelArray[],
                  TopGridData *MetaData, int FOFOnly);

/* From the problem bridge: reach the live hierarchy on a session handle. */
extern "C" LevelHierarchyEntry **EnzoModulesProblemLevelArray(void *h);
extern "C" TopGridData *EnzoModulesProblemMetaData(void *h);

/* Captured catalogue from the most recent find_halos call (single session,
 * one result set at a time -- the getters read this back into Python). */
namespace {
struct EMHalo {
  int    len;                 /* number of particles                  */
  double mass, mvir, rvir;    /* total / virial mass [Msun], r200 [kpc] */
  double vrms, spin;          /* velocity dispersion [km/s], spin     */
  double cm[3], vel[3], am[3];/* centre of mass, mean velocity, ang. mom. */
};
std::vector<EMHalo> em_halos;
int em_nsub = 0;              /* total SUBFIND subgroups (0 if not run) */
}

extern "C" {

/* Run the inline FOF halo finder (and SUBFIND if do_subfind) on the session's
 * particles.  linking_length / min_size override HaloFinderLinkingLength /
 * HaloFinderMinimumSize when > 0.  Returns the number of halos found (>= 0), or
 * -1 on failure.  Captures the catalogue for the halo_count / halo_property /
 * subhalo_count getters below.  No-op returning 0 if there are no particles. */
int enzomodules_session_find_halos(void *h, int do_subfind,
                                   double linking_length, int min_size)
{
  em_halos.clear();
  em_nsub = 0;

  LevelHierarchyEntry **LevelArray = EnzoModulesProblemLevelArray(h);
  TopGridData *MetaData = EnzoModulesProblemMetaData(h);

  /* No particles -> nothing to find.  Check up front, before touching the
   * finder: FOF_Initialize/the pipeline allocate the arrays deallocate_all_memory
   * frees, so bailing out after FOF_Initialize would free unallocated arrays. */
  long long npart = 0;
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (LevelHierarchyEntry *Temp = LevelArray[level]; Temp;
         Temp = Temp->NextGridThisLevel)
      npart += Temp->GridData->ReturnNumberOfParticles();
  if (npart == 0)
    return 0;

  /* Zero-initialize: the finder's struct is assumed value-initialized in
   * places (e.g. SUBFIND only ever += NSubGroupsAll), so a stack FOFData with
   * garbage counters would corrupt the catalogue. */
  FOFData AllVars;
  memset(&AllVars, 0, sizeof(FOFData));

  AllVars.LinkLength  = (linking_length > 0.0)
                      ? linking_length : HaloFinderLinkingLength;
  AllVars.GroupMinLen = (min_size > 0) ? min_size : HaloFinderMinimumSize;

  /* set_units enforces DesDensityNgb (32) <= GroupMinLen and would ENZO_FAIL
   * (throw) otherwise -- clamp so a small min_size can't abort the session. */
  if (AllVars.GroupMinLen < 32)
    AllVars.GroupMinLen = 32;

  set_units(AllVars);          /* units + SUBFIND neighbour counts + Theta   */
  AllVars.Grid = 256;          /* coarse-grid dimension (mirrors FOF.C)      */
  AllVars.MaxPlacement = 4;

  /* Copy enzo's grid particles into the finder's data structure. */
  FOF_Initialize(MetaData, LevelArray, AllVars, false);

  /* Serial friends-of-friends (the NumberOfProcessors == 1 path of FOF.C). */
  marking(AllVars);
  init_coarse_grid(AllVars);
  link_local_slab(AllVars);
  find_minids(AllVars);
  compile_group_catalogue(AllVars);

  if (do_subfind) {
    subfind(AllVars, MetaData->CycleNumber, MetaData->Time);
    em_nsub = AllVars.NSubGroupsAll;
  }

  /* Read each group's properties with the finder's own routines. */
  for (int gr = AllVars.NgroupsAll - 1; gr >= 0; gr--) {
    int head = AllVars.GroupDatAll[gr].Tag;
    int len  = AllVars.GroupDatAll[gr].Len;
    FOF_particle_data *Pbuf = new FOF_particle_data[len];
    get_particles(0, head, len, Pbuf, AllVars);

    float cm[3], cmv[3], AM[3], mtot, mstars, mvir, rvir, vrms, spin;
    get_properties(AllVars, Pbuf, len, false, cm, cmv, &mtot, &mstars,
                   &mvir, &rvir, AM, &vrms, &spin);

    EMHalo hh;
    hh.len = len;  hh.mass = mtot;  hh.mvir = mvir;  hh.rvir = rvir;
    hh.vrms = vrms;  hh.spin = spin;
    for (int d = 0; d < 3; d++) {
      hh.cm[d] = cm[d];  hh.vel[d] = cmv[d];  hh.am[d] = AM[d];
    }
    em_halos.push_back(hh);
    delete [] Pbuf;
  }

  int n = AllVars.NgroupsAll;

  /* Restore the particles to their grids.  FOF_Initialize MOVED them off the
   * hierarchy into AllVars.P (MoveParticlesFOF); FOF_Finalize copies them back
   * and redistributes them across the grids -- without this the next call (and
   * the rest of the session) would see empty/corrupt grids.  FOFOnly == FALSE
   * does the full RebuildHierarchy-style redistribution; it also frees D.P
   * (deallocate_all_memory deliberately leaves D.P for it). */
  FOF_Finalize(AllVars, LevelArray, MetaData, FALSE);
  deallocate_all_memory(AllVars);
  return n;
}

/* Number of halos captured by the most recent find_halos call. */
int enzomodules_session_halo_count(void)
{ return (int)em_halos.size(); }

/* Total SUBFIND subgroups from the most recent find_halos(do_subfind=1). */
int enzomodules_session_subhalo_count(void)
{ return em_nsub; }

/* Property `which` of halo `i` (0-based): 0 len, 1 mass, 2..4 centre of mass,
 * 5..7 mean velocity, 8 vrms, 9 spin, 10 virial mass, 11 virial radius, 12..14
 * angular momentum.  Returns 0 for out-of-range indices. */
double enzomodules_session_halo_property(int i, int which)
{
  if (i < 0 || i >= (int)em_halos.size()) return 0.0;
  const EMHalo &hh = em_halos[i];
  switch (which) {
    case 0:  return (double)hh.len;
    case 1:  return hh.mass;
    case 2:  return hh.cm[0];
    case 3:  return hh.cm[1];
    case 4:  return hh.cm[2];
    case 5:  return hh.vel[0];
    case 6:  return hh.vel[1];
    case 7:  return hh.vel[2];
    case 8:  return hh.vrms;
    case 9:  return hh.spin;
    case 10: return hh.mvir;
    case 11: return hh.rvir;
    case 12: return hh.am[0];
    case 13: return hh.am[1];
    case 14: return hh.am[2];
    default: return 0.0;
  }
}

} /* extern "C" */
