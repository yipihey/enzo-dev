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

int SetDefaultGlobalValues(TopGridData &MetaData);
int InitializeRateData(FLOAT Time);
FLOAT FindCrossSection(int type, float energy);

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

} /* extern "C" */
