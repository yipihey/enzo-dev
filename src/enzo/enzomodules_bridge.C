/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- implementation
/
/  Thin extern "C" shims forwarding to the legacy Fortran kernels.  See
/  enzomodules_bridge.h for the contract.  Nothing numerical lives here.
/
************************************************************************/

#include "enzomodules_bridge.h"

/* Fortran symbol name mangling.
 *
 * This mirrors FORTRAN_NAME in macros_and_parameters.h (default branch,
 * line ~136: NAME##_).  We do not include that header so the bridge can
 * also be built standalone (see EnzoModules/deps/build_pilot.sh) without
 * pulling in the full Enzo configuration.  When compiled as part of the
 * main `make lib` build the trailing-underscore convention still matches
 * the default; override EM_FORTRAN_NAME if a machine config differs. */
#ifndef EM_FORTRAN_NAME
#define EM_FORTRAN_NAME(NAME) NAME##_
#endif

/* Precision macros so the reported sizes track the build configuration.
 * Defaults match the pilot build (double baryons, 32-bit ints). */
#ifdef CONFIG_BFLOAT_4
#define EM_R_BYTES 4
#else
#define EM_R_BYTES 8
#endif

#ifdef LARGE_INTS
#define EM_I_BYTES 8
#else
#define EM_I_BYTES 4
#endif

extern "C" {

/* Legacy Fortran prototype (all arguments by reference). */
void EM_FORTRAN_NAME(twoshock)(
    double *dls, double *drs, double *pls, double *prs,
    double *uls, double *urs,
    int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2,
    double *dt, double *gamma, double *pmin, int *ipresfree,
    double *pbar, double *ubar, int *gravity, double *grslice,
    int *idual, double *eta1);

int enzomodules_baryon_precision_bytes(void) { return EM_R_BYTES; }
int enzomodules_int_precision_bytes(void)    { return EM_I_BYTES; }

void enzomodules_twoshock(
    double *dls, double *drs, double *pls, double *prs,
    double *uls, double *urs,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    double dt, double gamma, double pmin, int ipresfree,
    double *pbar, double *ubar, int gravity, double *grslice,
    int idual, double eta1)
{
  EM_FORTRAN_NAME(twoshock)(
      dls, drs, pls, prs, uls, urs,
      &idim, &jdim, &i1, &i2, &j1, &j2,
      &dt, &gamma, &pmin, &ipresfree,
      pbar, ubar, &gravity, grslice, &idual, &eta1);
}

} /* extern "C" */
