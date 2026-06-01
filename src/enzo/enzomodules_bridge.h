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

#ifdef __cplusplus
}
#endif

#endif /* ENZOMODULES_BRIDGE_H */
