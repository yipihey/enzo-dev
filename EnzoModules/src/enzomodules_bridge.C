/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- implementation
/
/  Thin extern "C" shims forwarding to the legacy Fortran kernels.  See
/  enzomodules_bridge.h for the contract.  Nothing numerical lives here.
/
************************************************************************/

#include "enzomodules_bridge.h"
#include <vector>
#include <cstring>
#include <cstdio>

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

/* When built standalone (the pilot library, see deps/build_pilot.sh) we do
 * not link Enzo's c_message.C / MPI stack, so provide minimal definitions of
 * the Fortran error/warning hooks that the kernels reference via the
 * ERROR_MESSAGE/WARNING_MESSAGE macros.  Guarded so the full `make lib` build
 * (which supplies the real fc_error/fc_warning) is unaffected. */
#ifdef ENZOMODULES_STANDALONE
extern "C" {
/* Fortran: subroutine f_error(sourcefile, linenumber); the trailing `long`
 * is the hidden CHARACTER*(*) length argument. */
void EM_FORTRAN_NAME(f_error)(char *msg, int *line, long len) {
  std::fprintf(stderr, "[EnzoModules] ENZO ERROR (%.*s:%d)\n",
               (int)len, msg ? msg : "", line ? *line : -1);
}
void EM_FORTRAN_NAME(f_warning)(char *msg, int *line, long len) {
  std::fprintf(stderr, "[EnzoModules] ENZO WARNING (%.*s:%d)\n",
               (int)len, msg ? msg : "", line ? *line : -1);
}
}
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

/* Legacy Fortran prototypes (all arguments by reference). */
void EM_FORTRAN_NAME(twoshock)(
    double *dls, double *drs, double *pls, double *prs,
    double *uls, double *urs,
    int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2,
    double *dt, double *gamma, double *pmin, int *ipresfree,
    double *pbar, double *ubar, int *gravity, double *grslice,
    int *idual, double *eta1);

void EM_FORTRAN_NAME(inteuler)(
    double *dslice, double *pslice, int *gravity, double *grslice,
    double *geslice, double *uslice, double *vslice, double *wslice,
    double *dxi, double *flatten,
    int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2,
    int *idual, double *eta1, double *eta2,
    int *isteep, int *iflatten, int *iconsrec, int *iposrec,
    double *dt, double *gamma, int *ipresfree,
    double *dls, double *drs, double *pls, double *prs,
    double *gels, double *gers, double *uls, double *urs,
    double *vls, double *vrs, double *wls, double *wrs,
    int *ncolor, double *colslice, double *colls, double *colrs);

void EM_FORTRAN_NAME(flux_twoshock)(
    double *dslice, double *eslice, double *geslice,
    double *uslice, double *vslice, double *wslice,
    double *dx, double *diffcoef,
    int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2,
    double *dt, double *gamma, int *idiff, int *idual, double *eta1,
    int *ifallback,
    double *dls, double *drs, double *pls, double *prs,
    double *gels, double *gers, double *uls, double *urs,
    double *vls, double *vrs, double *wls, double *wrs,
    double *pbar, double *ubar,
    double *df, double *ef, double *uf, double *vf, double *wf,
    double *gef, double *ges,
    int *ncolor, double *colslice, double *colls, double *colrs, double *colf);

void EM_FORTRAN_NAME(euler)(
    double *dslice, double *eslice, double *grslice, double *geslice,
    double *uslice, double *vslice, double *wslice,
    double *dx, double *diffcoef,
    int *idim, int *jdim, int *i1, int *i2, int *j1, int *j2,
    double *dt, double *gamma, int *idiff, int *gravity,
    int *idual, double *eta1, double *eta2,
    double *df, double *ef, double *uf, double *vf, double *wf,
    double *gef, double *ges,
    int *ncolor, double *colslice, double *colf, double *dfloor);

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

int enzomodules_ppm_sweep_1d(
    double *dslice, double *eslice,
    double *uslice, double *vslice, double *wslice,
    double *pslice,
    int idim, int i1, int i2,
    double dx, double dt, double gamma,
    double *df_out, double *ef_out, double *uf_out)
{
  /* This mirrors the body of grid::xEulerSweep (Grid_xEulerSweep.C) for the
     simplest faithful configuration.  Everything that is "off" in that
     configuration is still allocated and passed (the legacy kernels expect
     valid pointers); the corresponding feature flags are zero so the code
     paths are inert. */
  const int jdim = 1, j1 = 1, j2 = 1;
  const int n = idim * jdim;
  const int ie_p1 = i2 + 1;            /* twoshock solves one extra interface */

  /* Disabled features. */
  int    gravity = 0, idual = 0, ncolor = 0;
  int    isteep = 0, iflatten = 0, iconsrec = 0, iposrec = 0;
  int    ipresfree = 0, idiff = 0, ifallback = 0;
  double eta1 = 0.0, eta2 = 0.0;
  double pmin = 1e-20, dfloor = 1e-20;

  /* Scratch slabs (zero-initialised). */
  std::vector<double> grslice(n, 0.0), geslice(n, 0.0);
  std::vector<double> flatten(n, 0.0), diffcoef(n, 0.0);
  std::vector<double> dxa(n, dx);
  std::vector<double> dls(n, 0.0), drs(n, 0.0), pls(n, 0.0), prs(n, 0.0);
  std::vector<double> uls(n, 0.0), urs(n, 0.0), vls(n, 0.0), vrs(n, 0.0);
  std::vector<double> wls(n, 0.0), wrs(n, 0.0), gels(n, 0.0), gers(n, 0.0);
  std::vector<double> pbar(n, 0.0), ubar(n, 0.0);
  std::vector<double> df(n, 0.0), ef(n, 0.0), uf(n, 0.0), vf(n, 0.0);
  std::vector<double> wf(n, 0.0), gef(n, 0.0), ges(n, 0.0);
  /* Colour scratch: length-n dummies (ncolor == 0 -> loops are inert). */
  std::vector<double> colslice(n, 0.0), colls(n, 0.0), colrs(n, 0.0), colf(n, 0.0);

  /* 1. PPM reconstruction -> left/right interface states. */
  EM_FORTRAN_NAME(inteuler)(
      dslice, pslice, &gravity, grslice.data(), geslice.data(),
      uslice, vslice, wslice, dxa.data(), flatten.data(),
      (int*)&idim, (int*)&jdim, &i1, &i2, (int*)&j1, (int*)&j2,
      &idual, &eta1, &eta2, &isteep, &iflatten, &iconsrec, &iposrec,
      &dt, &gamma, &ipresfree,
      dls.data(), drs.data(), pls.data(), prs.data(),
      gels.data(), gers.data(), uls.data(), urs.data(),
      vls.data(), vrs.data(), wls.data(), wrs.data(),
      &ncolor, colslice.data(), colls.data(), colrs.data());

  /* 2. Lagrangian Riemann problem at each interface. */
  EM_FORTRAN_NAME(twoshock)(
      dls.data(), drs.data(), pls.data(), prs.data(),
      uls.data(), urs.data(),
      (int*)&idim, (int*)&jdim, &i1, (int*)&ie_p1, (int*)&j1, (int*)&j2,
      &dt, &gamma, &pmin, &ipresfree,
      pbar.data(), ubar.data(), &gravity, grslice.data(), &idual, &eta1);

  /* 3. Eulerian fluxes from the resolved states. */
  EM_FORTRAN_NAME(flux_twoshock)(
      dslice, eslice, geslice.data(), uslice, vslice, wslice,
      dxa.data(), diffcoef.data(),
      (int*)&idim, (int*)&jdim, &i1, &i2, (int*)&j1, (int*)&j2,
      &dt, &gamma, &idiff, &idual, &eta1, &ifallback,
      dls.data(), drs.data(), pls.data(), prs.data(),
      gels.data(), gers.data(), uls.data(), urs.data(),
      vls.data(), vrs.data(), wls.data(), wrs.data(),
      pbar.data(), ubar.data(),
      df.data(), ef.data(), uf.data(), vf.data(), wf.data(),
      gef.data(), ges.data(),
      &ncolor, colslice.data(), colls.data(), colrs.data(), colf.data());

  /* 4. Conservative update of the zone-centred quantities (in place). */
  EM_FORTRAN_NAME(euler)(
      dslice, eslice, grslice.data(), geslice.data(),
      uslice, vslice, wslice, dxa.data(), diffcoef.data(),
      (int*)&idim, (int*)&jdim, &i1, &i2, (int*)&j1, (int*)&j2,
      &dt, &gamma, &idiff, &gravity, &idual, &eta1, &eta2,
      df.data(), ef.data(), uf.data(), vf.data(), wf.data(),
      gef.data(), ges.data(),
      &ncolor, colslice.data(), colf.data(), &dfloor);

  if (df_out) std::memcpy(df_out, df.data(), n * sizeof(double));
  if (ef_out) std::memcpy(ef_out, ef.data(), n * sizeof(double));
  if (uf_out) std::memcpy(uf_out, uf.data(), n * sizeof(double));
  return 0;
}

} /* extern "C" */
