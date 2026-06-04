/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE
/
/  PURPOSE:
/    Stable, extern "C" surface that exposes selected legacy Enzo compute
/    kernels to external callers (in particular EnzoModules.jl via ccall).
/    Each entry point is a thin shim that forwards to the existing Fortran
/    kernel; no behaviour is reimplemented here.
/
/    This is modelled on the libyt integration (ExposeHierarchyToLibyt.C):
/    a narrow C ABI over the in-tree implementation, compiled into the
/    shared library so downstream code can call it without linking C++.
/
/    ABI/precision contract: the caller MUST match the precision this
/    library was built with.  Query it at runtime via the
/    enzomodules_*_precision_bytes() functions rather than assuming.
/
************************************************************************/

#ifndef ENZOMODULES_BRIDGE_H
#define ENZOMODULES_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Size in bytes of the baryon floating-point type (R_PREC) this library
   was compiled with: 4 for CONFIG_BFLOAT_4, 8 for CONFIG_BFLOAT_8. */
int enzomodules_baryon_precision_bytes(void);

/* Size in bytes of the Fortran integer type (INTG_PREC) this library was
   compiled with: 4 for SMALL_INTS, 8 for LARGE_INTS. */
int enzomodules_int_precision_bytes(void);

/* ------------------------------------------------------------------ *
 *  Hydro: two-shock approximate Riemann solver  (twoshock.F)
 *
 *  Slicewise solver over a 2D (idim x jdim) column-major slab.  Arrays
 *  hold the reconstructed left/right interface states; pbar/ubar receive
 *  the resolved interface pressure and normal velocity.  Indices i1..i2,
 *  j1..j2 are 1-based (Fortran) and inclusive.
 *
 *  Scalars are passed by value here and forwarded by reference to the
 *  Fortran routine.  Note pls/prs may be overwritten by the kernel when
 *  ipresfree != 0 (matching the legacy in/out semantics).
 * ------------------------------------------------------------------ */
void enzomodules_twoshock(
    double *dls, double *drs,        /* left/right density            */
    double *pls, double *prs,        /* left/right pressure  (in/out) */
    double *uls, double *urs,        /* left/right normal velocity    */
    int idim, int jdim,              /* declared slab dimensions      */
    int i1, int i2, int j1, int j2,  /* 1-based inclusive loop bounds */
    double dt, double gamma, double pmin, int ipresfree,
    double *pbar, double *ubar,      /* resolved interface state (out)*/
    int gravity, double *grslice, int idual, double eta1);

/* ------------------------------------------------------------------ *
 *  Hydro: one directional PPM hydro update of a 1D slice
 *
 *  Composite kernel that reproduces the numerical core of the Enzo
 *  Eulerian PPM solver's directional sweep (src/enzo/Grid_xEulerSweep.C):
 *
 *      inteuler  ->  twoshock  ->  flux_twoshock  ->  euler
 *
 *  for the simplest faithful configuration: no gravity, no dual energy,
 *  no colour fields, no artificial diffusion / slope flattening, PPM
 *  reconstruction + two-shock Riemann solver.  It is the worked example
 *  for wrapping a *composite* solver (not just a leaf kernel) behind the
 *  bridge.
 *
 *  The slice arrays are 1-based Fortran slabs of length idim (jdim = 1).
 *  Indices i1..i2 are the active (non-ghost) cells; PPM needs >= 3 ghost
 *  cells on each side, so the caller must size idim = (i2-i1+1) + 2*NGHOST
 *  with NGHOST >= 3 and fill/refresh ghosts (boundary conditions) itself.
 *
 *  Inputs  : dslice (density), eslice (total specific energy),
 *            uslice/vslice/wslice (velocities), pslice (pressure,
 *            precomputed by the caller -- Enzo computes it once per step).
 *  Outputs : dslice, eslice, uslice, vslice, wslice are updated in place
 *            with the new state after a timestep dt.  Optionally df/ef/uf
 *            (density/energy/normal-momentum fluxes) are returned if the
 *            pointers are non-NULL.
 *
 *  dx is the (uniform) cell width.  Returns 0 on success.
 * ------------------------------------------------------------------ */
int enzomodules_ppm_sweep_1d(
    double *dslice, double *eslice,
    double *uslice, double *vslice, double *wslice,
    double *pslice,
    int idim, int i1, int i2,
    double dx, double dt, double gamma,
    double *df, double *ef, double *uf);   /* flux outputs, may be NULL */

/* ------------------------------------------------------------------ *
 *  Individual PPM-stage kernels (component-level certification).
 *
 *  One shim per Fortran leaf kernel, with every production feature
 *  (dual energy via geslice/eta1/eta2, gravity via grslice, colour via
 *  ncolor/colslice, flattening/diffusion) exposed as a parameter.  Arrays
 *  are caller-allocated column-major (idim x jdim) slabs; in/out semantics
 *  follow the underlying routine.  See enzomodules_bridge.C for details.
 * ------------------------------------------------------------------ */
void enzomodules_pgas2d(
    double *dslice, double *eslice, double *pslice,
    double *uslice, double *vslice, double *wslice,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    double gamma, double pmin);

void enzomodules_pgas2d_dual(
    double *dslice, double *eslice, double *geslice, double *pslice,
    double *uslice, double *vslice, double *wslice,
    double eta1, double eta2,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    double gamma, double pmin);

void enzomodules_calcdiss(
    double *dslice, double *eslice, double *uslice, double *v, double *w,
    double *pslice, double *dx, double *dy, double *dz,
    int idim, int jdim, int kdim, int i1, int i2, int j1, int j2,
    int k, int nzz, int idir, int dimx, int dimy, int dimz,
    double dt, double gamma, int idiff, int iflatten,
    double *diffcoef, double *flatten);

void enzomodules_inteuler(
    double *dslice, double *pslice, int gravity, double *grslice,
    double *geslice, double *uslice, double *vslice, double *wslice,
    double *dxi, double *flatten,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    int idual, double eta1, double eta2,
    int isteep, int iflatten, int iconsrec, int iposrec,
    double dt, double gamma, int ipresfree,
    double *dls, double *drs, double *pls, double *prs,
    double *gels, double *gers, double *uls, double *urs,
    double *vls, double *vrs, double *wls, double *wrs,
    int ncolor, double *colslice, double *colls, double *colrs);

void enzomodules_flux_twoshock(
    double *dslice, double *eslice, double *geslice,
    double *uslice, double *vslice, double *wslice,
    double *dx, double *diffcoef,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    double dt, double gamma, int idiff, int idual, double eta1,
    int ifallback,
    double *dls, double *drs, double *pls, double *prs,
    double *gels, double *gers, double *uls, double *urs,
    double *vls, double *vrs, double *wls, double *wrs,
    double *pbar, double *ubar,
    double *df, double *ef, double *uf, double *vf, double *wf,
    double *gef, double *ges,
    int ncolor, double *colslice, double *colls, double *colrs, double *colf);

void enzomodules_euler(
    double *dslice, double *eslice, double *grslice, double *geslice,
    double *uslice, double *vslice, double *wslice,
    double *dx, double *diffcoef,
    int idim, int jdim, int i1, int i2, int j1, int j2,
    double dt, double gamma, int idiff, int gravity,
    int idual, double eta1, double eta2,
    double *df, double *ef, double *uf, double *vf, double *wf,
    double *gef, double *ges,
    int ncolor, double *colslice, double *colf, double dfloor);

/* Full production directional sweep (dual energy + gravity + colour +
   flattening/diffusion as parameters).  geslice/grslice/colslice may be
   NULL when the corresponding feature is off.  Returns 0 on success. */
int enzomodules_ppm_sweep_1d_full(
    double *dslice, double *eslice, double *geslice,
    double *uslice, double *vslice, double *wslice, double *pslice,
    int idim, int i1, int i2, double dx, double dt, double gamma,
    int gravity, double *grslice,
    int idual, double eta1, double eta2,
    int isteep, int iflatten, int iconsrec, int iposrec,
    int idiff, int ipresfree, int ifallback,
    double pmin, double dfloor,
    int ncolor, double *colslice,
    double *df, double *ef, double *uf);

#ifdef __cplusplus
}
#endif

#endif /* ENZOMODULES_BRIDGE_H */
