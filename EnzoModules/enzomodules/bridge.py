"""ctypes loader + raw bindings for the EnzoModules pilot shared library.

Locates and loads ``libenzomodules_pilot.so`` (built by
``deps/build_pilot.sh``) and exposes the C ABI declared in
``src/enzo/enzomodules_bridge.h``.

Library resolution order:
  1. ``$ENZOMODULES_LIB`` if set,
  2. ``deps/libenzomodules_pilot.so`` next to this package.

If the library is absent, :func:`available` returns ``False`` and the higher
layers skip live calls, so a checkout without the C build can still parse
fixtures and exercise the diff utilities.
"""
from __future__ import annotations

import ctypes
import os
from typing import Optional, Sequence, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_LIB = os.path.join(_HERE, "..", "deps", "libenzomodules_pilot.so")

_lib: Optional[ctypes.CDLL] = None


def libpath() -> str:
    """Absolute path to the pilot shared library (may not exist)."""
    env = os.environ.get("ENZOMODULES_LIB", "")
    return os.path.abspath(env) if env else os.path.abspath(_DEFAULT_LIB)


def available() -> bool:
    """True if the pilot shared library exists on disk."""
    return os.path.isfile(libpath())


def _load() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        if not available():
            raise RuntimeError(
                f"EnzoModules pilot library not found at {libpath()}.\n"
                "Build it with EnzoModules/deps/build_pilot.sh, or set "
                "ENZOMODULES_LIB.")
        lib = ctypes.CDLL(libpath())
        lib.enzomodules_baryon_precision_bytes.restype = ctypes.c_int
        lib.enzomodules_baryon_precision_bytes.argtypes = []
        lib.enzomodules_int_precision_bytes.restype = ctypes.c_int
        lib.enzomodules_int_precision_bytes.argtypes = []
        d = ctypes.POINTER(ctypes.c_double)
        lib.enzomodules_twoshock.restype = None
        lib.enzomodules_twoshock.argtypes = [
            d, d, d, d, d, d,                                   # dls drs pls prs uls urs
            ctypes.c_int, ctypes.c_int,                         # idim jdim
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,  # i1 i2 j1 j2
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int,  # dt gamma pmin ipresfree
            d, d,                                               # pbar ubar
            ctypes.c_int, d, ctypes.c_int, ctypes.c_double,     # gravity grslice idual eta1
        ]
        lib.enzomodules_ppm_sweep_1d.restype = ctypes.c_int
        lib.enzomodules_ppm_sweep_1d.argtypes = [
            d, d, d, d, d, d,                                   # dslice eslice uslice vslice wslice pslice
            ctypes.c_int, ctypes.c_int, ctypes.c_int,          # idim i1 i2
            ctypes.c_double, ctypes.c_double, ctypes.c_double,  # dx dt gamma
            d, d, d,                                            # df ef uf (may be NULL)
        ]
        _lib = lib
    return _lib


def check_precision() -> Tuple[int, int]:
    """Verify the library precision matches the bindings (double, int32).

    Returns ``(baryon_bytes, int_bytes)``; raises on mismatch -- precisely the
    class of silent bug this framework exists to catch.
    """
    lib = _load()
    rb = lib.enzomodules_baryon_precision_bytes()
    ib = lib.enzomodules_int_precision_bytes()
    if rb != ctypes.sizeof(ctypes.c_double):
        raise RuntimeError(
            f"EnzoModules built with {rb}-byte baryons; bindings assume "
            f"{ctypes.sizeof(ctypes.c_double)}.")
    if ib != ctypes.sizeof(ctypes.c_int):
        raise RuntimeError(
            f"EnzoModules built with {ib}-byte Fortran ints; bindings assume "
            f"{ctypes.sizeof(ctypes.c_int)}.")
    return rb, ib


def _carr(vals: Sequence[float]) -> ctypes.Array:
    return (ctypes.c_double * len(vals))(*[float(v) for v in vals])


def twoshock_raw(dls, drs, pls, prs, uls, urs,
                 idim, jdim, i1, i2, j1, j2,
                 dt, gamma, pmin, ipresfree,
                 gravity, grslice, idual, eta1):
    """Direct binding to ``enzomodules_twoshock``.

    Sequences are length ``idim*jdim`` (column-major).  Returns
    ``(pbar, ubar)`` as plain Python lists.
    """
    lib = _load()
    n = idim * jdim
    c_dls, c_drs = _carr(dls), _carr(drs)
    c_pls, c_prs = _carr(pls), _carr(prs)
    c_uls, c_urs = _carr(uls), _carr(urs)
    c_grslice = _carr(grslice)
    c_pbar = (ctypes.c_double * n)()
    c_ubar = (ctypes.c_double * n)()
    lib.enzomodules_twoshock(
        c_dls, c_drs, c_pls, c_prs, c_uls, c_urs,
        int(idim), int(jdim), int(i1), int(i2), int(j1), int(j2),
        float(dt), float(gamma), float(pmin), int(ipresfree),
        c_pbar, c_ubar,
        int(gravity), c_grslice, int(idual), float(eta1),
    )
    return list(c_pbar), list(c_ubar)


# --------------------------------------------------------------------------
# hydro_rk (Runge-Kutta MUSCL) line solvers, incl. Dedner-cleaning MHD.
#
# These live in a *separate* shared library (libenzomodules_hydrork.so) which
# links against the full Enzo shared library, because the C++ solvers depend
# on Enzo's headers and global state (unlike the standalone Fortran kernels).
# --------------------------------------------------------------------------

_HYDRORK_DEFAULT = os.path.join(_HERE, "..", "deps", "libenzomodules_hydrork.so")
_hydrork = None

#: Riemann solver codes (Enzo typedefs.h enum).
HLL, LLF, HLLC, HLLD = 1, 3, 4, 6


def hydrork_libpath():
    env = os.environ.get("ENZOMODULES_HYDRORK_LIB", "")
    return os.path.abspath(env) if env else os.path.abspath(_HYDRORK_DEFAULT)


def hydrork_available():
    return os.path.isfile(hydrork_libpath())


def _load_hydrork():
    global _hydrork
    if _hydrork is None:
        if not hydrork_available():
            raise RuntimeError(
                f"hydro_rk library not found at {hydrork_libpath()}.\n"
                "Build it with EnzoModules/deps/build_hydro_rk.sh (needs the "
                "full Enzo shared library), or set ENZOMODULES_HYDRORK_LIB.")
        lib = ctypes.CDLL(hydrork_libpath())
        d = ctypes.POINTER(ctypes.c_double)
        common = [ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_int,
                  ctypes.c_double, ctypes.c_double]
        lib.enzomodules_hydro_rk_line.restype = ctypes.c_int
        lib.enzomodules_hydro_rk_line.argtypes = common + [
            d, ctypes.c_int, ctypes.c_int, ctypes.c_int, d]
        lib.enzomodules_mhd_rk_line.restype = ctypes.c_int
        lib.enzomodules_mhd_rk_line.argtypes = common + [
            ctypes.c_double, d, ctypes.c_int, ctypes.c_int, ctypes.c_int, d]
        _hydrork = lib
    return _hydrork


def _line_flux(func, prim_rows, active_size, extra, gamma, theta_limiter,
               nghost, small_rho, small_p):
    neq = len(prim_rows)
    ncells = len(prim_rows[0])
    flat = [v for row in prim_rows for v in row]
    c_prim = (ctypes.c_double * len(flat))(*[float(v) for v in flat])
    c_flux = (ctypes.c_double * (neq * (active_size + 1)))()
    args = [HLL, gamma, theta_limiter, nghost, small_rho, small_p]  # solver filled by caller
    rc = func(*extra, c_prim, neq, ncells, active_size, c_flux)
    if rc != 0:
        raise RuntimeError(f"{func.__name__} returned {rc}")
    return [[c_flux[f * (active_size + 1) + i] for i in range(active_size + 1)]
            for f in range(neq)]


def hydro_rk_line(solver, prim_rows, gamma=1.4, theta_limiter=1.5, nghost=3,
                  small_rho=1e-20, small_p=1e-20):
    """Reconstruction + Riemann flux for one hydro line (NEQ_HYDRO=5).

    ``prim_rows`` is 5 rows ``[rho, eint, vx, vy, vz]`` each of length
    ``ncells = active + 2*nghost``.  Returns 5 flux rows of length ``active+1``.
    """
    lib = _load_hydrork()
    active = len(prim_rows[0]) - 2 * nghost
    extra = [solver, gamma, theta_limiter, nghost, small_rho, small_p]
    return _line_flux(lib.enzomodules_hydro_rk_line, prim_rows, active, extra,
                      gamma, theta_limiter, nghost, small_rho, small_p)


def mhd_rk_line(solver, prim_rows, c_h, gamma=1.4, theta_limiter=1.5, nghost=3,
                small_rho=1e-20, small_p=1e-20):
    """Reconstruction + Riemann flux for one Dedner-MHD line (NEQ_MHD=9).

    ``prim_rows`` is 9 rows ``[rho, eint, vx, vy, vz, Bx, By, Bz, Phi]``.
    ``c_h`` is the Dedner divergence-cleaning wave speed.  Returns 9 flux rows.
    """
    lib = _load_hydrork()
    active = len(prim_rows[0]) - 2 * nghost
    extra = [solver, gamma, theta_limiter, nghost, small_rho, small_p, c_h]
    return _line_flux(lib.enzomodules_mhd_rk_line, prim_rows, active, extra,
                      gamma, theta_limiter, nghost, small_rho, small_p)


# --------------------------------------------------------------------------
# Grid-method solvers (ZEUS, and later radiation transfer, gravity).
#
# These are grid:: *methods*; the bridge builds a minimal grid fixture and
# calls the method.  Like the hydro_rk library, libenzomodules_grid.so links
# against the full Enzo shared library.
# --------------------------------------------------------------------------

_GRID_DEFAULT = os.path.join(_HERE, "..", "deps", "libenzomodules_grid.so")
_gridlib = None


def grid_libpath():
    env = os.environ.get("ENZOMODULES_GRID_LIB", "")
    return os.path.abspath(env) if env else os.path.abspath(_GRID_DEFAULT)


def grid_available():
    return os.path.isfile(grid_libpath())


def _load_gridlib():
    global _gridlib
    if _gridlib is None:
        if not grid_available():
            raise RuntimeError(
                f"grid solver library not found at {grid_libpath()}.\n"
                "Build it with EnzoModules/deps/build_grid.sh (needs the full "
                "Enzo shared library), or set ENZOMODULES_GRID_LIB.")
        lib = ctypes.CDLL(grid_libpath())
        d = ctypes.POINTER(ctypes.c_double)
        lib.enzomodules_zeus_sweep_1d.restype = ctypes.c_int
        lib.enzomodules_zeus_sweep_1d.argtypes = [
            d, d, d, ctypes.c_int, ctypes.c_int,
            ctypes.c_double, ctypes.c_double, ctypes.c_double]
        _gridlib = lib
    return _gridlib


def cic_deposit(particles, dims, left=(0.0, 0.0, 0.0), right=(1.0, 1.0, 1.0)):
    """CIC-deposit particles onto a 3D grid (legacy grid::DepositParticlePositions).

    ``particles`` is a list of ``(x, y, z, mass)``.  ``dims`` is the active
    ``(nx, ny, nz)``.  Returns ``(field, cellvol)`` where ``field`` is the
    deposited gravitating-mass field (a flat list, including the gravity
    buffer) and ``cellvol`` is its cell volume.  The deposit is linear and
    conservative; the absolute mass scale follows Enzo's particle-mass unit
    convention, so tests use ratios / invariants.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_cic_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_cic_deposit.restype = ctypes.c_int
        lib.enzomodules_cic_deposit.argtypes = [
            ctypes.c_int, ip, d, d, d, d, d, d, ctypes.c_int,
            d, ctypes.c_int, ip, d]
        lib._cic_set = True
    n = len(particles)
    c_dims = (ctypes.c_int * 3)(int(dims[0]), int(dims[1]), int(dims[2]))
    c_left = (ctypes.c_double * 3)(*[float(v) for v in left])
    c_right = (ctypes.c_double * 3)(*[float(v) for v in right])
    px = (ctypes.c_double * n)(*[float(p[0]) for p in particles])
    py = (ctypes.c_double * n)(*[float(p[1]) for p in particles])
    pz = (ctypes.c_double * n)(*[float(p[2]) for p in particles])
    mass = (ctypes.c_double * n)(*[float(p[3]) for p in particles])
    cap = (max(dims) + 20) ** 3
    out = (ctypes.c_double * cap)()
    osize = ctypes.c_int(0)
    ocv = ctypes.c_double(0.0)
    rc = lib.enzomodules_cic_deposit(3, c_dims, c_left, c_right, px, py, pz,
                                     mass, n, out, cap, ctypes.byref(osize),
                                     ctypes.byref(ocv))
    if rc != 0:
        raise RuntimeError(f"enzomodules_cic_deposit returned {rc} "
                           f"(size {osize.value} > capacity {cap}?)")
    return [out[k] for k in range(osize.value)], ocv.value


# Canonical species ordering by MultiSpecies level (matches the C bridge).
SPECIES_6 = ["De", "HI", "HII", "HeI", "HeII", "HeIII"]
SPECIES_9 = SPECIES_6 + ["HM", "H2I", "H2II"]
SPECIES_12 = SPECIES_9 + ["DI", "DII", "HDI"]
SPECIES_BY_LEVEL = {1: SPECIES_6, 2: SPECIES_9, 3: SPECIES_12}


def chemistry_step(rho, e_tot, species, dt, multispecies=None,
                   density_units=1.673e-24, length_units=3.086e21,
                   time_units=3.156e13):
    """Advance one primordial-chemistry + radiative-cooling step (legacy
    grid::SolveRateAndCoolEquations, the non-Grackle MultiSpecies network).

    Supports all three MultiSpecies levels:
      1 ->  6 species (De, HI, HII, HeI, HeII, HeIII)
      2 ->  9 species (+ HM, H2I, H2II)
      3 -> 12 species (+ DI, DII, HDI)

    ``species`` is a dict {name: list} of code-unit mass densities; ``rho`` and
    ``e_tot`` are lists.  ``multispecies`` is inferred from the species present
    if not given.  Returns ``(e_tot, species_dict)`` updated.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_chem_set"):
        d = ctypes.POINTER(ctypes.c_double)
        lib.enzomodules_chemistry_step.restype = ctypes.c_int
        lib.enzomodules_chemistry_step.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
            ctypes.c_double, ctypes.c_double, d, d, d]
        lib.enzomodules_num_species.restype = ctypes.c_int
        lib.enzomodules_num_species.argtypes = [ctypes.c_int]
        lib._chem_set = True

    if multispecies is None:
        multispecies = 3 if "DI" in species else 2 if "H2I" in species else 1
    names = SPECIES_BY_LEVEL[multispecies]
    n = len(rho)

    flat = []
    for name in names:
        flat.extend(float(v) for v in species[name])
    c_flat = (ctypes.c_double * (len(names) * n))(*flat)
    c_rho = (ctypes.c_double * n)(*[float(v) for v in rho])
    c_e = (ctypes.c_double * n)(*[float(v) for v in e_tot])

    rc = lib.enzomodules_chemistry_step(
        int(multispecies), n, float(dt), float(density_units),
        float(length_units), float(time_units), c_rho, c_e, c_flat)
    if rc != 0:
        raise RuntimeError(f"enzomodules_chemistry_step returned {rc}")

    out = {name: [c_flat[s * n + i] for i in range(n)]
           for s, name in enumerate(names)}
    return [c_e[i] for i in range(n)], out


def poisson_solve(rhs, dims):
    """Solve the discrete Poisson equation L(phi) = rhs with Enzo's multigrid
    solver (the engine behind grid::SolveForPotential).

    ``rhs`` is a flat list of length ``prod(dims)``; ``dims`` is ``(nx, ny, nz)``
    with trailing 1s for lower rank (use 2^k+1 per active axis for clean
    coarsening).  Returns ``(solution, norm, mean)`` where ``norm/mean`` is the
    converged residual diagnostic.  The operator carries an internal cell-size
    normalization (see SolveForPotential's ``Constant``), so manufactured-
    solution checks compare shape.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_poisson_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_poisson_solve.restype = ctypes.c_int
        lib.enzomodules_poisson_solve.argtypes = [ctypes.c_int, ip, d, d, d, d]
        lib._poisson_set = True
    rank = sum(1 for v in dims if v > 1) or 1
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    size = len(rhs)
    c_rhs = (ctypes.c_double * size)(*[float(v) for v in rhs])
    c_sol = (ctypes.c_double * size)()
    norm = ctypes.c_double(0.0)
    mean = ctypes.c_double(0.0)
    rc = lib.enzomodules_poisson_solve(rank, c_dims, c_rhs, c_sol,
                                       ctypes.byref(norm), ctypes.byref(mean))
    if rc != 0:
        raise RuntimeError(f"enzomodules_poisson_solve returned {rc}")
    return [c_sol[i] for i in range(size)], norm.value, mean.value


def hi_cross_section(energy):
    """Enzo's HI photo-ionization cross-section [cm^2] at ``energy`` [eV]
    (Verner et al. 1996 fits)."""
    lib = _load_gridlib()
    if not hasattr(lib, "_xsec_set"):
        lib.enzomodules_hi_cross_section.restype = ctypes.c_double
        lib.enzomodules_hi_cross_section.argtypes = [ctypes.c_double]
        lib._xsec_set = True
    return lib.enzomodules_hi_cross_section(float(energy))


def raytrace_uniform(hi_density, energy=14.0, photons=1e50, path_fraction=0.3,
                     nx=48, density_units=1.673e-24, length_units=3.086e21,
                     time_units=3.156e13):
    """Trace one photon package through a uniform-HI grid (legacy
    grid::WalkPhotonPackage) and return ``(photons_final, radius, kph_sum)``.

    With ``density_units = m_H`` the HI number density is ``hi_density`` cm^-3.
    The ray travels ``path_fraction`` box lengths; verify Beer-Lambert with
    ``photons_final = photons * exp(-hi_density * sigma * radius*length_units)``
    using :func:`hi_cross_section`.  A negative ``photons_final`` means the ray
    was optically thick enough to be fully absorbed.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_ray_set"):
        d = ctypes.POINTER(ctypes.c_double)
        lib.enzomodules_raytrace_uniform.restype = ctypes.c_int
        lib.enzomodules_raytrace_uniform.argtypes = [
            ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            d, d, d]
        lib._ray_set = True
    pf = (ctypes.c_double * 1)()
    rf = (ctypes.c_double * 1)()
    ks = (ctypes.c_double * 1)()
    rc = lib.enzomodules_raytrace_uniform(
        int(nx), float(density_units), float(length_units), float(time_units),
        float(hi_density), float(energy), float(photons), float(path_fraction),
        pf, rf, ks)
    if rc != 0:
        raise RuntimeError(f"enzomodules_raytrace_uniform returned {rc}")
    return pf[0], rf[0], ks[0]


# AMR refinement criteria (CellFlaggingMethod codes) exercisable here.
FLAG_SLOPE = 1
FLAG_MASS = 2              # baryon mass / overdensity
FLAG_SHOCKS = 3
FLAG_SHEAR = 9
FLAG_SECOND_DERIVATIVE = 15


def flag_cells(rank, dims, density, method=FLAG_SLOPE, threshold=0.3,
               energy=None, u=None, v=None, w=None, dx=1.0):
    """AMR cell flagging via the real grid::SetFlaggingField dispatch, for a
    chosen refinement criterion (``method``, a CellFlaggingMethod code):
    1=slope, 2=baryon mass/overdensity, 3=shocks, 9=shear,
    15=second derivative.  ``density`` is required; ``energy``/``u``/``v``/``w``
    default to a quiescent state (needed by shock/shear criteria).  Returns
    ``(flagging, count)``.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_flag_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_flag_cells.restype = ctypes.c_int
        lib.enzomodules_flag_cells.argtypes = [
            ctypes.c_int, ip, ctypes.c_double, ctypes.c_int, ctypes.c_double,
            d, d, d, d, d, ip, ip]
        lib._flag_set = True
    size = len(density)
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    zero = [0.0] * size
    energy = zero if energy is None else energy
    u = zero if u is None else u
    v = zero if v is None else v
    w = zero if w is None else w

    def arr(a):
        return (ctypes.c_double * size)(*[float(x) for x in a])

    c_d, c_e = arr(density), arr(energy)
    c_u, c_v, c_w = arr(u), arr(v), arr(w)
    c_flag = (ctypes.c_int * size)()
    c_count = ctypes.c_int(0)
    rc = lib.enzomodules_flag_cells(int(rank), c_dims, float(dx), int(method),
                                    float(threshold), c_d, c_e, c_u, c_v, c_w,
                                    c_flag, ctypes.byref(c_count))
    if rc != 0:
        raise RuntimeError(f"enzomodules_flag_cells returned {rc}")
    return [c_flag[i] for i in range(size)], c_count.value


def flag_region(dims, region_left, region_right, dx=1.0):
    """AMR must-refine-region flagging (CellFlaggingMethod 12, 3D only) -- flag
    every cell inside the box [region_left, region_right] (domain coords).
    Returns ``(flagging, count)``."""
    lib = _load_gridlib()
    if not hasattr(lib, "_region_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_flag_region.restype = ctypes.c_int
        lib.enzomodules_flag_region.argtypes = [
            ctypes.c_int, ip, ctypes.c_double, d, d, ip, ip]
        lib._region_set = True
    c_dims = (ctypes.c_int * 3)(int(dims[0]), int(dims[1]), int(dims[2]))
    rl = (ctypes.c_double * 3)(*[float(x) for x in region_left])
    rr = (ctypes.c_double * 3)(*[float(x) for x in region_right])
    size = dims[0] * dims[1] * dims[2]
    c_flag = (ctypes.c_int * size)()
    c_count = ctypes.c_int(0)
    rc = lib.enzomodules_flag_region(3, c_dims, float(dx), rl, rr, c_flag,
                                     ctypes.byref(c_count))
    if rc != 0:
        raise RuntimeError(f"enzomodules_flag_region returned {rc}")
    return [c_flag[i] for i in range(size)], c_count.value


def flag_jeans(rank, dims, density, energy, safety=4.0,
               density_units=1.673e-24, length_units=3.086e21,
               time_units=3.156e13, dx=1.0):
    """AMR Jeans-length flagging (CellFlaggingMethod 6) -- refine where the
    Jeans length is under-resolved (dense/cold gas).  Returns
    ``(flagging, count)``."""
    lib = _load_gridlib()
    if not hasattr(lib, "_jeans_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_flag_jeans.restype = ctypes.c_int
        lib.enzomodules_flag_jeans.argtypes = [
            ctypes.c_int, ip, ctypes.c_double, d, d, ctypes.c_double,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ip, ip]
        lib._jeans_set = True
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    size = len(density)
    c_d = (ctypes.c_double * size)(*[float(x) for x in density])
    c_e = (ctypes.c_double * size)(*[float(x) for x in energy])
    c_flag = (ctypes.c_int * size)()
    c_count = ctypes.c_int(0)
    rc = lib.enzomodules_flag_jeans(int(rank), c_dims, float(dx), c_d, c_e,
                                    float(safety), float(density_units),
                                    float(length_units), float(time_units),
                                    c_flag, ctypes.byref(c_count))
    if rc != 0:
        raise RuntimeError(f"enzomodules_flag_jeans returned {rc}")
    return [c_flag[i] for i in range(size)], c_count.value


def cluster(rank, dims, flagging, dx=1.0, max_sub=200):
    """Cluster a flagging field into rectangular subgrids (Berger-Rigoutsos
    via ProtoSubgrid + IdentifyNewSubgridsBySignature).

    ``flagging`` is the int FlaggingField (length prod(dims)).  Returns a list
    of subgrid boxes, each ``[(xlo,xhi),(ylo,yhi),(zlo,zhi)]`` in active-cell
    indices -- one per generated child grid.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_cluster_set"):
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_cluster.restype = ctypes.c_int
        lib.enzomodules_cluster.argtypes = [
            ctypes.c_int, ip, ctypes.c_double, ip, ip, ip, ctypes.c_int]
        lib._cluster_set = True
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    c_flag = (ctypes.c_int * len(flagging))(*[int(x) for x in flagging])
    c_nsub = ctypes.c_int(0)
    c_boxes = (ctypes.c_int * (6 * max_sub))()
    rc = lib.enzomodules_cluster(int(rank), c_dims, float(dx), c_flag,
                                 ctypes.byref(c_nsub), c_boxes, int(max_sub))
    if rc != 0:
        raise RuntimeError(f"enzomodules_cluster returned {rc}")
    n = min(c_nsub.value, max_sub)
    return [[(c_boxes[s * 6 + 2 * d], c_boxes[s * 6 + 2 * d + 1])
             for d in range(3)] for s in range(n)]


def ppm_hydro_step(rank, dims, d, e, u, v, w, dx, dt, gamma):
    """One full multi-dimensional PPM hydro step (legacy
    grid::SolveHydroEquations -> x/y/z Euler sweeps).

    Builds a minimal Enzo grid of the given ``rank`` and active+ghost ``dims``
    (domain edges ``[0, dims[k]*dx]``), sets the Density/TotalEnergy/Velocity
    fields, takes one in-place PPM_DirectEuler step over ``dt``, and reads the
    fields back.  ``d`` (density), ``e`` (specific *total* energy = internal +
    0.5 v^2), ``u``/``v``/``w`` (Velocity1/2/3) are flat row-major lists of
    length ``prod(dims)`` (INCLUDING ghost zones).  Returns updated
    ``(d, e, u, v, w)`` as new lists.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_ppm_step_set"):
        d_p = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_ppm_hydro_step.restype = ctypes.c_int
        lib.enzomodules_ppm_hydro_step.argtypes = [
            ctypes.c_int, ip, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            d_p, d_p, d_p, d_p, d_p]
        lib._ppm_step_set = True
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    size = len(d)
    cd = (ctypes.c_double * size)(*[float(x) for x in d])
    ce = (ctypes.c_double * size)(*[float(x) for x in e])
    cu = (ctypes.c_double * size)(*[float(x) for x in u])
    cv = (ctypes.c_double * size)(*[float(x) for x in v])
    cw = (ctypes.c_double * size)(*[float(x) for x in w])
    rc = lib.enzomodules_ppm_hydro_step(int(rank), c_dims, float(dx),
                                        float(dt), float(gamma),
                                        cd, ce, cu, cv, cw)
    if rc != 0:
        raise RuntimeError(f"enzomodules_ppm_hydro_step returned {rc}")
    return list(cd), list(ce), list(cu), list(cv), list(cw)


def mhdct_step(dims, d, e, u, v, w, bx, by, bz, dx, dt, gamma, nsteps=1):
    """``nsteps`` constrained-transport (CT) MHD steps (grid::SolveMHD_Li via
    grid::SolveHydroEquations, HydroMethod=MHD_Li, UseMHDCT=1).

    Lays down both the cell-centered Bfield1/2/3 and the staggered face-centered
    MagneticField[], takes ``nsteps`` in-place CT steps on the same grid (so the
    face-centered field carries over), and reads the fields back.  Arrays are
    flat row-major of length prod(dims) incl. ghost zones.  Returns
    ``(d, e, u, v, w, bx, by, bz, maxdivb_before, maxdivb_after)`` -- the CT
    invariant is that both divergences stay ~machine precision.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_mhdct_step_set"):
        d_p = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_mhdct_step.restype = ctypes.c_int
        lib.enzomodules_mhdct_step.argtypes = [
            ip, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int,
            d_p, d_p, d_p, d_p, d_p, d_p, d_p, d_p, d_p, d_p]
        lib._mhdct_step_set = True
    c_dims = (ctypes.c_int * 3)(int(dims[0]),
                                int(dims[1]) if len(dims) > 1 else 1,
                                int(dims[2]) if len(dims) > 2 else 1)
    size = len(d)
    cd = (ctypes.c_double * size)(*[float(x) for x in d])
    ce = (ctypes.c_double * size)(*[float(x) for x in e])
    cu = (ctypes.c_double * size)(*[float(x) for x in u])
    cv = (ctypes.c_double * size)(*[float(x) for x in v])
    cw = (ctypes.c_double * size)(*[float(x) for x in w])
    cbx = (ctypes.c_double * size)(*[float(x) for x in bx])
    cby = (ctypes.c_double * size)(*[float(x) for x in by])
    cbz = (ctypes.c_double * size)(*[float(x) for x in bz])
    before = ctypes.c_double(0.0)
    after = ctypes.c_double(0.0)
    rc = lib.enzomodules_mhdct_step(c_dims, float(dx), float(dt), float(gamma),
                                    int(nsteps),
                                    cd, ce, cu, cv, cw, cbx, cby, cbz,
                                    ctypes.byref(before), ctypes.byref(after))
    if rc != 0:
        raise RuntimeError(f"enzomodules_mhdct_step returned {rc}")
    return (list(cd), list(ce), list(cu), list(cv), list(cw),
            list(cbx), list(cby), list(cbz), before.value, after.value)


def raytrace_twogrid(hi_density, energy=14.0, photons=1e50, path_fraction=0.6,
                     nx=32, density_units=1.673e-24, length_units=3.086e21,
                     time_units=3.156e13):
    """Multi-grid (AMR) photon transport: trace a ray across two grids tiling
    the domain in x (legacy SubgridMarker -> FindPhotonNewGrid handoff).

    Returns ``(photons_final, radius, grids_visited)``.  ``grids_visited == 2``
    means the ray crossed the grid boundary; attenuation should still obey
    Beer-Lambert for the total path (verifying continuity across grids).
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_ray2_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_raytrace_twogrid.restype = ctypes.c_int
        lib.enzomodules_raytrace_twogrid.argtypes = [
            ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
            d, d, ip]
        lib._ray2_set = True
    pf = (ctypes.c_double * 1)()
    rf = (ctypes.c_double * 1)()
    gv = (ctypes.c_int * 1)()
    rc = lib.enzomodules_raytrace_twogrid(
        int(nx), float(density_units), float(length_units), float(time_units),
        float(hi_density), float(energy), float(photons), float(path_fraction),
        pf, rf, gv)
    if rc != 0:
        raise RuntimeError(f"enzomodules_raytrace_twogrid returned {rc}")
    return pf[0], rf[0], gv[0]


def rt_identify(kph_field, nghost=3, dx=0.05):
    """Build a full grid carrying the radiative-transfer rate fields and run
    Enzo's field identification (grid::IdentifyRadiativeTransferFields).

    Proves the full-grid fixture supports radiation-transport data structures.
    ``kph_field`` is the kphHI (photo-ionization rate) field; returns
    ``(kphHI_out, kphHINum, gammaNum)`` — the round-tripped field and the
    located field indices.
    """
    lib = _load_gridlib()
    if not hasattr(lib, "_rt_set"):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_rt_identify.restype = ctypes.c_int
        lib.enzomodules_rt_identify.argtypes = [
            ctypes.c_int, ctypes.c_int, ctypes.c_double, d, d, ip, ip]
        lib._rt_set = True
    idim = len(kph_field)
    cin = (ctypes.c_double * idim)(*[float(v) for v in kph_field])
    cout = (ctypes.c_double * idim)()
    kn = ctypes.c_int(0)
    gn = ctypes.c_int(0)
    rc = lib.enzomodules_rt_identify(idim, nghost, float(dx), cin, cout,
                                     ctypes.byref(kn), ctypes.byref(gn))
    if rc != 0:
        raise RuntimeError(f"enzomodules_rt_identify returned {rc}")
    return [cout[i] for i in range(idim)], kn.value, gn.value


def zeus_sweep_1d(d_arr, e_arr, u_arr, nghost, dx, dt, gamma):
    """One ZEUS hydro update of a 1D slice (legacy grid::ZeusSolver).

    ``d_arr`` (density), ``e_arr`` (specific *internal* energy), ``u_arr``
    (x-velocity) are length ``idim = active + 2*nghost``.  Returns updated
    ``(d, e, u)`` as new lists.
    """
    lib = _load_gridlib()
    idim = len(d_arr)
    cd = (ctypes.c_double * idim)(*[float(v) for v in d_arr])
    ce = (ctypes.c_double * idim)(*[float(v) for v in e_arr])
    cu = (ctypes.c_double * idim)(*[float(v) for v in u_arr])
    rc = lib.enzomodules_zeus_sweep_1d(cd, ce, cu, idim, nghost,
                                       float(dx), float(dt), float(gamma))
    if rc != 0:
        raise RuntimeError(f"enzomodules_zeus_sweep_1d returned {rc}")
    return list(cd), list(ce), list(cu)


def ppm_sweep_1d(dslice, eslice, uslice, vslice, wslice, pslice,
                 i1, i2, dx, dt, gamma, want_fluxes=False):
    """Direct binding to ``enzomodules_ppm_sweep_1d``.

    Runs one PPM hydro update of a 1D slice (inteuler -> twoshock ->
    flux_twoshock -> euler).  Sequences are length ``idim``; ``i1``/``i2`` are
    the 1-based inclusive active-cell range (>= 3 ghost cells each side).
    Returns the updated ``(dslice, eslice, uslice, vslice, wslice)`` as new
    lists, and -- if ``want_fluxes`` -- a ``(df, ef, uf)`` tuple, else None.
    """
    lib = _load()
    idim = len(dslice)
    c_d, c_e = _carr(dslice), _carr(eslice)
    c_u, c_v, c_w = _carr(uslice), _carr(vslice), _carr(wslice)
    c_p = _carr(pslice)
    if want_fluxes:
        c_df = (ctypes.c_double * idim)()
        c_ef = (ctypes.c_double * idim)()
        c_uf = (ctypes.c_double * idim)()
    else:
        c_df = c_ef = c_uf = None
    rc = lib.enzomodules_ppm_sweep_1d(
        c_d, c_e, c_u, c_v, c_w, c_p,
        int(idim), int(i1), int(i2),
        float(dx), float(dt), float(gamma),
        c_df, c_ef, c_uf,
    )
    if rc != 0:
        raise RuntimeError(f"enzomodules_ppm_sweep_1d returned {rc}")
    state = (list(c_d), list(c_e), list(c_u), list(c_v), list(c_w))
    fluxes = (list(c_df), list(c_ef), list(c_uf)) if want_fluxes else None
    return state, fluxes
