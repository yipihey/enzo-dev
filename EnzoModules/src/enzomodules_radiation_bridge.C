/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- radiative transfer (photon ray-tracing)
/
/  PURPOSE:
/    Fire a single photon package through a uniform-density grid and report the
/    surviving photon count and path length, so callers can verify
/    Beer-Lambert attenuation  N(L) = N0 * exp(-n_HI * sigma * L)  end-to-end
/    through the legacy ray-tracer (grid::WalkPhotonPackage).  Also exposes
/    Enzo's HI photo-ionization cross-section (FindCrossSection) so the test
/    can compute the expected optical depth independently.
/
/    Compiled with Enzo headers, linked against libenzo.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "RadiativeTransferHealpixRoutines64.h"

int SetDefaultGlobalValues(TopGridData &MetaData);
int InitializeRateData(FLOAT Time);
FLOAT FindCrossSection(int type, float energy);

/* Set the radiative-transfer globals for a uniform-medium ray trace. */
static void em_set_rt_globals(double density_units, double length_units,
                              double time_units)
{
  RadiativeTransfer                         = 1;
  RadiativeTransferHydrogenOnly             = 1;
  MultiSpecies                              = 1;
  RadiationFieldType                        = 0;
  ComovingCoordinates                       = 0;
  RadiationPressure                         = 0;
  RadiativeTransferSourceClustering         = 0;
  RadiativeTransferSplitPhotonRadius        = 1e-30;
  RadiativeTransferAdaptiveTimestep         = 0;
  RadiativeTransferPropagationSpeedFraction = 1.0;
  RadiativeTransferPhotonEscapeRadius       = 0.0;
  RadiativeTransferOpticallyThinH2          = 0;
  NumberOfGhostZones                        = 3;
  GlobalDensityUnits = density_units;
  GlobalLengthUnits  = length_units;
  GlobalTimeUnits    = time_units;
  GlobalMassUnits    = density_units * length_units * length_units * length_units;
  for (int d = 0; d < MAX_DIMENSION; d++) {
    DomainLeftEdge[d]  = 0.0;
    DomainRightEdge[d] = 1.0;
  }
  PhotonTime = 0.0;
}

/* HEALPix pixel (at the given level) whose direction is closest to +x. */
static int em_best_px_for_plus_x(int level)
{
  int64_t nside = (int64_t)(1 << level);
  int64_t npix = 12 * nside * nside;
  int best = 0;
  double bestdot = -2.0, vec[3];
  for (int64_t p = 0; p < npix; p++) {
    pix2vec_nest64(nside, p, vec);
    if (vec[0] > bestdot) { bestdot = vec[0]; best = (int)p; }
  }
  return best;
}

extern "C" {

/* Enzo's HI photo-ionization cross-section [cm^2] at `energy` [eV]
 * (Verner et al. 1996 fits) -- the same one the ray-tracer uses. */
double enzomodules_hi_cross_section(double energy)
{
  return (double) FindCrossSection(0, (float)energy);
}

/* Trace one photon package of `photons` photons at `energy` eV through a
 * uniform grid (nx^3) with HI number density `hi_density` (code units; with
 * density_units = m_H this is n_HI in cm^-3).  The ray travels `path_fraction`
 * box lengths (periodic wrap keeps it in the uniform medium).  Returns the
 * surviving photons, path length (code units), and summed kphHI.  Returns 0
 * on success. */
int enzomodules_raytrace_uniform(
    int nx, double density_units, double length_units, double time_units,
    double hi_density, double energy, double photons, double path_fraction,
    double *photons_final, double *radius_final, double *kph_sum)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);

  RadiativeTransfer                         = 1;
  RadiativeTransferHydrogenOnly             = 1;
  MultiSpecies                              = 1;
  RadiationFieldType                        = 0;
  ComovingCoordinates                       = 0;
  RadiationPressure                         = 0;
  RadiativeTransferSourceClustering         = 0;
  RadiativeTransferSplitPhotonRadius        = 1e-30;   /* disable splitting */
  RadiativeTransferAdaptiveTimestep         = 0;
  RadiativeTransferPropagationSpeedFraction = 1.0;
  RadiativeTransferPhotonEscapeRadius       = 0.0;
  RadiativeTransferOpticallyThinH2          = 0;
  NumberOfGhostZones                        = 3;

  GlobalDensityUnits = density_units;
  GlobalLengthUnits  = length_units;
  GlobalTimeUnits    = time_units;
  GlobalMassUnits    = density_units * length_units * length_units * length_units;

  for (int d = 0; d < MAX_DIMENSION; d++) {
    DomainLeftEdge[d]  = 0.0;
    DomainRightEdge[d] = 1.0;
  }
  PhotonTime = 0.0;

  InitializeRateData(0.0);

  double velocity_units = length_units / time_units;
  double lightspeed     = clight / velocity_units;          /* code units */
  double dtphoton       = path_fraction / lightspeed;       /* travel path_fraction box-lengths */
  dtPhoton              = dtphoton;

  grid g;
  int   dims[3]  = { nx, nx, nx };
  FLOAT left[3]  = { (FLOAT)0.0, (FLOAT)0.0, (FLOAT)0.0 };
  FLOAT right[3] = { (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0 };
  int   ftypes[13] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3,
                       ElectronDensity, HIDensity, HIIDensity,
                       HeIDensity, HeIIDensity, HeIIIDensity, kphHI, PhotoGamma };
  int size = g.EnzoModulesSetupGrid(3, dims, left, right, 13, ftypes, dtphoton);

  std::vector<double> ones(size, 1.0), tiny(size, 1e-20), zero(size, 0.0);
  std::vector<double> hi(size, hi_density), egas(size, 1e-3);
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Density),         ones.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(TotalEnergy),     egas.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity1),       zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity2),       zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(Velocity3),       zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HIDensity),       hi.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HIIDensity),      tiny.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(ElectronDensity), tiny.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIDensity),      tiny.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIIDensity),     tiny.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(HeIIIDensity),    tiny.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(kphHI),           zero.data());
  g.EnzoModulesSetField(g.EnzoModulesFieldIndex(PhotoGamma),      zero.data());

  return g.EnzoModulesRaytrace(energy, photons, dtphoton, lightspeed,
                               /*ipix*/ 0, /*level*/ 0,
                               photons_final, radius_final, kph_sum);
}

/* Multi-grid (AMR) photon transport: trace a ray through two grids tiling the
 * domain in x (gA = [0,0.5], gB = [0.5,1]), both uniform HI.  The ray is fired
 * nearly along +x so it crosses from gA into gB.  Returns the surviving
 * photons, total path length, and the number of grids the ray visited (2 means
 * it successfully crossed the grid boundary).  Returns 0 on success. */
int enzomodules_raytrace_twogrid(
    int nx, double density_units, double length_units, double time_units,
    double hi_density, double energy, double photons, double path_fraction,
    double *photons_final, double *radius_final, int *grids_visited)
{
  TopGridData MetaData;
  SetDefaultGlobalValues(MetaData);
  em_set_rt_globals(density_units, length_units, time_units);
  InitializeRateData(0.0);

  double velocity_units = length_units / time_units;
  double lightspeed     = clight / velocity_units;
  double dtphoton       = path_fraction / lightspeed;
  dtPhoton              = dtphoton;

  const int hlevel = 3;                       /* Nside=8: a near-+x direction */
  const int ipix   = em_best_px_for_plus_x(hlevel);
  const int ng     = NumberOfGhostZones;

  int ftypes[13] = { Density, TotalEnergy, Velocity1, Velocity2, Velocity3,
                     ElectronDensity, HIDensity, HIIDensity,
                     HeIDensity, HeIIDensity, HeIIIDensity, kphHI, PhotoGamma };

  /* gA: x in [0,0.5]; gB: x in [0.5,1]; both full in y,z.  Equal cell size. */
  int dimsA[3] = { nx / 2 + 2 * ng, nx + 2 * ng, nx + 2 * ng };
  grid gA, gB;
  FLOAT lA[3] = { 0.0, 0.0, 0.0 }, rA[3] = { 0.5, 1.0, 1.0 };
  FLOAT lB[3] = { 0.5, 0.0, 0.0 }, rB[3] = { 1.0, 1.0, 1.0 };
  int sizeA = gA.EnzoModulesSetupGrid(3, dimsA, lA, rA, 13, ftypes, dtphoton);
  int sizeB = gB.EnzoModulesSetupGrid(3, dimsA, lB, rB, 13, ftypes, dtphoton);

  std::vector<double> zA(sizeA, 0.0), tA(sizeA, 1e-20), hA(sizeA, hi_density), eA(sizeA, 1e-3), oA(sizeA, 1.0);
  std::vector<double> zB(sizeB, 0.0), tB(sizeB, 1e-20), hB(sizeB, hi_density), eB(sizeB, 1e-3), oB(sizeB, 1.0);
  grid *gs[2] = { &gA, &gB };
  std::vector<double> *z[2] = { &zA, &zB }, *t[2] = { &tA, &tB },
                      *h[2] = { &hA, &hB }, *e[2] = { &eA, &eB }, *o[2] = { &oA, &oB };
  for (int k = 0; k < 2; k++) {
    grid *g = gs[k];
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(Density),         o[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(TotalEnergy),     e[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(Velocity1),       z[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(Velocity2),       z[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(Velocity3),       z[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(HIDensity),       h[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(HIIDensity),      t[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(ElectronDensity), t[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(HeIDensity),      t[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(HeIIDensity),     t[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(HeIIIDensity),    t[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(kphHI),           z[k]->data());
    g->EnzoModulesSetField(g->EnzoModulesFieldIndex(PhotoGamma),      z[k]->data());
  }

  return gA.EnzoModulesRaytraceTwoGrid(&gB, energy, photons, dtphoton,
                                       lightspeed, ipix, hlevel,
                                       photons_final, radius_final,
                                       grids_visited);
}

} /* extern "C" */
