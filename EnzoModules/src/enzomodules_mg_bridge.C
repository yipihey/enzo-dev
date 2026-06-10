/***********************************************************************
/
/  ENZOMODULES C-ABI BRIDGE -- multigrid Poisson kernels
/
/  PURPOSE:
/    Expose Enzo's multigrid Poisson Fortran kernels (mg_relax, mg_calc_defect,
/    mg_restrict, mg_prolong) and the acceleration gradient (comp_accel) as
/    flat-array C-ABI oracles, so the KernelAbstractions port (lib/PoissonKernels)
/    can be certified bit-for-bit against the live legacy kernels.  Nothing
/    numerical is reimplemented here -- these are thin marshalling wrappers around
/    the Fortran routines, which live in libenzo.
/
/    PRECISION: this bridge is linked against the p8_b8 Enzo library, where
/    R_PREC = real*8 (double) and INTG_PREC = integer*4 (int).  The multigrid
/    Fortran routines therefore take double* data and int* dims -- declared so
/    below.  (The float* prototypes in MultigridSolver.C are for a b4 build only.)
/
/    LAYOUT: arrays are Fortran column-major, flat length prod(dims); Julia
/    Array{Float64,3} shares that layout exactly, so the pointers pass straight
/    through with no transpose.
/
************************************************************************/

#include "macros_and_parameters.h"

/* Multigrid Fortran kernels (defined in libenzo; b8 => double data, int dims). */
extern "C" {
void FORTRAN_NAME(mg_relax)(double *solution, double *rhs, int *ndim,
                            int *dim1, int *dim2, int *dim3);
void FORTRAN_NAME(mg_calc_defect)(double *solution, double *rhs, double *defect,
                            int *ndim, int *dim1, int *dim2, int *dim3,
                            double *norm);
void FORTRAN_NAME(mg_restrict)(double *source, double *dest, int *ndim,
                            int *sdim1, int *sdim2, int *sdim3,
                            int *ddim1, int *ddim2, int *ddim3);
void FORTRAN_NAME(mg_prolong)(double *source, double *dest, int *ndim,
                            int *sdim1, int *sdim2, int *sdim3,
                            int *ddim1, int *ddim2, int *ddim3);
void FORTRAN_NAME(comp_accel)(double *source, double *dest1, double *dest2,
                            double *dest3, int *ndim, int *iflag,
                            int *sdim1, int *sdim2, int *sdim3,
                            int *ddim1, int *ddim2, int *ddim3,
                            int *start1, int *start2, int *start3,
                            double *delx, double *dely, double *delz);
}

extern "C" {

/* One Gauss-Seidel relaxation (in place on `solution`).  Returns 0. */
int enzomodules_mg_relax(int ndim, int dim1, int dim2, int dim3,
                         double *solution, double *rhs)
{
  FORTRAN_NAME(mg_relax)(solution, rhs, &ndim, &dim1, &dim2, &dim3);
  return 0;
}

/* Negative residual into `defect` + its L2 norm into *out_norm.  Returns 0. */
int enzomodules_mg_calc_defect(int ndim, int dim1, int dim2, int dim3,
                               double *solution, double *rhs,
                               double *defect, double *out_norm)
{
  double norm = 0.0;
  FORTRAN_NAME(mg_calc_defect)(solution, rhs, defect, &ndim,
                               &dim1, &dim2, &dim3, &norm);
  *out_norm = norm;
  return 0;
}

/* Restrict fine `source` (sdims) onto coarse `dest` (ddims).  Returns 0. */
int enzomodules_mg_restrict(int ndim, int sdim1, int sdim2, int sdim3,
                            int ddim1, int ddim2, int ddim3,
                            double *source, double *dest)
{
  FORTRAN_NAME(mg_restrict)(source, dest, &ndim,
                            &sdim1, &sdim2, &sdim3, &ddim1, &ddim2, &ddim3);
  return 0;
}

/* Prolong coarse `source` (sdims) onto fine `dest` (ddims).  Returns 0. */
int enzomodules_mg_prolong(int ndim, int sdim1, int sdim2, int sdim3,
                           int ddim1, int ddim2, int ddim3,
                           double *source, double *dest)
{
  FORTRAN_NAME(mg_prolong)(source, dest, &ndim,
                           &sdim1, &sdim2, &sdim3, &ddim1, &ddim2, &ddim3);
  return 0;
}

/* Difference potential `source` (sdims) into the three acceleration fields
 * (ddims), with dest->source offset (start1..3), staggering iflag, cell sizes
 * del{x,y,z}.  Returns 0. */
int enzomodules_comp_accel(int ndim, int iflag,
                           int sdim1, int sdim2, int sdim3,
                           int ddim1, int ddim2, int ddim3,
                           int start1, int start2, int start3,
                           double delx, double dely, double delz,
                           double *source, double *dest1, double *dest2,
                           double *dest3)
{
  FORTRAN_NAME(comp_accel)(source, dest1, dest2, dest3, &ndim, &iflag,
                           &sdim1, &sdim2, &sdim3, &ddim1, &ddim2, &ddim3,
                           &start1, &start2, &start3, &delx, &dely, &delz);
  return 0;
}

} /* extern "C" */
