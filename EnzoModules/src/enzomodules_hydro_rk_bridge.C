/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- hydro_rk (Runge-Kutta MUSCL) solvers
/
/  PURPOSE:
/    Expose the hydro_rk 1D line Riemann solvers (reconstruction + flux) to
/    external callers, for both hydrodynamics and Dedner-cleaning MHD.  These
/    are free functions in Enzo (HydroLine / MHDLine dispatch to HLL_PLM,
/    HLLC_PLM, LLF_PLM / HLLD_PLM_MHD, ...), so the bridge sets the handful of
/    globals they read and calls them -- nothing numerical is reimplemented.
/
/    Unlike the Fortran-kernel bridge (enzomodules_bridge.C, built standalone),
/    this translation unit is compiled against the full Enzo headers and linked
/    against the Enzo shared library, which provides the globals' definitions.
/
/    Primitive layout (rows):
/      hydro (NEQ_HYDRO=5): [ rho, eint, vx, vy, vz ]
/      MHD   (NEQ_MHD  =9): [ rho, eint, vx, vy, vz, Bx, By, Bz, Phi ]
/    where eint is *internal* specific energy.  Arrays are row-major
/    [field*ncells + cell]; the line has `active_size` active cells plus
/    `nghost` ghost cells on each side (so ncells = active_size + 2*nghost).
/    Fluxes are returned as [field*(active_size+1) + interface].
/
************************************************************************/

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* hydro_rk line dispatchers (declared here; defined in libenzo). */
int HydroLine(float **Prim, float **priml, float **primr,
              float **species, float **colors, float **FluxLine, int ActiveSize,
              float dtdx, char direc, int ij, int ik, int fallback);
int MHDLine(float **Prim, float **priml, float **primr,
            float **species, float **colors, float **FluxLine, int ActiveSize,
            float dtdx, char direc, int jj, int kk, int fallback);

/* Build a float** view (nrows pointers into a freshly copied flat buffer). */
static float **alloc_rows(int nrows, int ncols, const double *flat) {
  float **rows = new float *[nrows];
  for (int f = 0; f < nrows; f++) {
    rows[f] = new float[ncols];
    if (flat) for (int i = 0; i < ncols; i++) rows[f][i] = (float)flat[f * ncols + i];
    else      for (int i = 0; i < ncols; i++) rows[f][i] = 0.0;
  }
  return rows;
}
static void free_rows(float **rows, int nrows) {
  for (int f = 0; f < nrows; f++) delete[] rows[f];
  delete[] rows;
}

/* Common globals shared by hydro + MHD line solvers. */
static void set_common_globals(int neq, double gamma, double theta_limiter,
                               int nghost, double small_rho, double small_p) {
  Gamma              = gamma;
  EOSType            = 0;            /* ideal gas */
  Theta_Limiter      = theta_limiter;
  NumberOfGhostZones = nghost;
  DualEnergyFormalism = 0;
  NSpecies = 0;
  NColor   = 0;
  SmallRho = small_rho;
  SmallP   = small_p;
  SmallT   = 0.0;
  SmallEint = 0.0;
  Mu       = 0.6;
  /* primitive/flux index conventions (SetDefaultGlobalValues defaults) */
  iden = 0; ieint = 0; ivx = 1; ivy = 2; ivz = 3; ietot = 4;
  iBx = 5; iBy = 6; iBz = 7; iPhi = 8;
  iD = 0; iEint = 0; iS1 = 1; iS2 = 2; iS3 = 3; iEtot = 4;
}

extern "C" {

/* Run one hydro line: reconstruction + Riemann flux.
 * riemann_solver: HLL(1) | LLF(3) | HLLC(4).  Returns 0 on success. */
int enzomodules_hydro_rk_line(
    int riemann_solver, double gamma, double theta_limiter, int nghost,
    double small_rho, double small_p,
    const double *prim_flat, int neq, int ncells, int active_size,
    double *flux_flat)
{
  set_common_globals(neq, gamma, theta_limiter, nghost, small_rho, small_p);
  NEQ_HYDRO = neq;
  RiemannSolver = riemann_solver;
  ReconstructionMethod = PLM;

  float **Prim  = alloc_rows(neq, ncells, prim_flat);
  float **priml = alloc_rows(neq, active_size + 1, 0);
  float **primr = alloc_rows(neq, active_size + 1, 0);
  float **flux  = alloc_rows(neq, active_size + 1, 0);
  float **species = alloc_rows(1, active_size + 1, 0);
  float **colors  = alloc_rows(1, active_size + 1, 0);

  int rc = HydroLine(Prim, priml, primr, species, colors, flux, active_size,
                     /*dtdx*/ 0.0, 'x', 0, 0, /*fallback*/ 0);

  if (rc != FAIL)
    for (int f = 0; f < neq; f++)
      for (int i = 0; i < active_size + 1; i++)
        flux_flat[f * (active_size + 1) + i] = (double)flux[f][i];

  free_rows(Prim, neq); free_rows(priml, neq); free_rows(primr, neq);
  free_rows(flux, neq); free_rows(species, 1); free_rows(colors, 1);
  return (rc == FAIL) ? 1 : 0;
}

/* Run one MHD line with Dedner divergence cleaning.
 * riemann_solver: HLLD(6) | HLL(1) | LLF(3).  c_h is the Dedner wave speed.
 * Returns 0 on success. */
int enzomodules_mhd_rk_line(
    int riemann_solver, double gamma, double theta_limiter, int nghost,
    double small_rho, double small_p, double c_h,
    const double *prim_flat, int neq, int ncells, int active_size,
    double *flux_flat)
{
  set_common_globals(neq, gamma, theta_limiter, nghost, small_rho, small_p);
  NEQ_MHD = neq;
  RiemannSolver = riemann_solver;
  ReconstructionMethod = PLM;
  C_h = c_h;

  float **Prim  = alloc_rows(neq, ncells, prim_flat);
  float **priml = alloc_rows(neq, active_size + 1, 0);
  float **primr = alloc_rows(neq, active_size + 1, 0);
  float **flux  = alloc_rows(neq, active_size + 1, 0);
  float **species = alloc_rows(1, active_size + 1, 0);
  float **colors  = alloc_rows(1, active_size + 1, 0);

  int rc = MHDLine(Prim, priml, primr, species, colors, flux, active_size,
                   /*dtdx*/ 0.0, 'x', 0, 0, /*fallback*/ 0);

  if (rc != FAIL)
    for (int f = 0; f < neq; f++)
      for (int i = 0; i < active_size + 1; i++)
        flux_flat[f * (active_size + 1) + i] = (double)flux[f][i];

  free_rows(Prim, neq); free_rows(priml, neq); free_rows(primr, neq);
  free_rows(flux, neq); free_rows(species, 1); free_rows(colors, 1);
  return (rc == FAIL) ? 1 : 0;
}

} /* extern "C" */
