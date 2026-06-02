"""Generate initial conditions for any Enzo problem type from its parameter
file, by driving Enzo's own ``InitializeNew`` (which dispatches on
``ProblemType`` to every problem generator) and exposing the resulting grid
hierarchy.

    from enzomodules.problems import Problem
    with Problem("Toro-1-ShockTube.enzo") as p:
        g = p.grid(0)
        rho = g.field("Density")          # flat list, includes ghost zones

Requires the grid solver library (deps/build_grid.sh); ``available()`` reports
whether it is built.
"""
from __future__ import annotations

import ctypes
import os
import tempfile
from typing import List, Optional

from . import bridge

# A practical subset of Enzo FieldType enum values (typedefs.h) for naming.
FIELD_NAMES = {
    0: "Density", 1: "TotalEnergy", 2: "InternalEnergy", 3: "Pressure",
    4: "Velocity1", 5: "Velocity2", 6: "Velocity3", 7: "ElectronDensity",
    8: "HIDensity", 9: "HIIDensity", 10: "HeIDensity", 11: "HeIIDensity",
    12: "HeIIIDensity", 13: "HMDensity", 14: "H2IDensity", 15: "H2IIDensity",
    19: "Metallicity", 20: "GravPotential",
    23: "kphHI", 24: "PhotoGamma",
    50: "Bfield1", 51: "Bfield2", 52: "Bfield3", 53: "PhiField",
}
NAME_TO_FIELD = {v: k for k, v in FIELD_NAMES.items()}


import contextlib


@contextlib.contextmanager
def _suppress_fd_output():
    """Redirect C-level stdout/stderr (fd 1/2) to /dev/null -- EvolveHierarchy
    prints a per-cycle log we don't want in test output."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = (os.dup(1), os.dup(2))
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(saved[0])
        os.close(saved[1])
        os.close(devnull)


def available() -> bool:
    return bridge.grid_available()


def _lib():
    lib = bridge._load_gridlib()
    if not getattr(lib, "_problem_set", False):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_init_problem.restype = ctypes.c_void_p
        lib.enzomodules_init_problem.argtypes = [ctypes.c_char_p]
        lib.enzomodules_evolve_problem.restype = ctypes.c_void_p
        lib.enzomodules_evolve_problem.argtypes = [
            ctypes.c_char_p, ctypes.c_double, ctypes.c_int]
        for fn in ("num_grids", "grid_rank", "num_fields", "grid_size",
                   "problemtype", "num_particles"):
            f = getattr(lib, "enzomodules_problem_" + fn)
            f.restype = ctypes.c_int
            f.argtypes = ([ctypes.c_void_p] if fn in ("num_grids", "problemtype")
                          else [ctypes.c_void_p, ctypes.c_int])
        lib.enzomodules_problem_grid_dims.argtypes = [ctypes.c_void_p, ctypes.c_int, ip]
        lib.enzomodules_problem_field_types.argtypes = [ctypes.c_void_p, ctypes.c_int, ip]
        lib.enzomodules_problem_get_field.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, d]
        lib.enzomodules_problem_get_particle_pos.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, d]
        lib.enzomodules_free_problem.argtypes = [ctypes.c_void_p]
        # Persistent-session stepping API.
        lib.enzomodules_session_init.restype = ctypes.c_void_p
        lib.enzomodules_session_init.argtypes = [ctypes.c_char_p]
        # Data output / restart.
        lib.enzomodules_session_write_output.restype = ctypes.c_int
        lib.enzomodules_session_write_output.argtypes = [
            ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        lib.enzomodules_session_from_output.restype = ctypes.c_void_p
        lib.enzomodules_session_from_output.argtypes = [ctypes.c_char_p]
        for fn in ("time", "stop_time", "compute_dt"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = ctypes.c_double
        lib.enzomodules_session_time.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_stop_time.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_cycle.restype = ctypes.c_int
        lib.enzomodules_session_cycle.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_compute_dt.argtypes = [ctypes.c_void_p, ctypes.c_int]
        for fn in ("set_boundary", "solve_hydro", "solve_cooling",
                   "star_particles", "rebuild", "gravity", "evolve_photons",
                   "create_fluxes", "update_from_finer", "finalize_fluxes",
                   "num_grids_on_level"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = ctypes.c_int
            f.argtypes = [ctypes.c_void_p, ctypes.c_int]
        for fn in ("set_dt", "advance_time", "update_particles",
                   "copy_baryon_to_old", "clear_boundary_fluxes"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = None
        lib.enzomodules_session_set_dt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        # EvolvePhotons with a source-mode flag (0 = parameter-file sources,
        # 1 = convert star particles into radiation sources).
        lib.enzomodules_session_evolve_photons_ex.restype = ctypes.c_int
        lib.enzomodules_session_evolve_photons_ex.argtypes = [
            ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
        # Inline FOF halo finder + SUBFIND.
        lib.enzomodules_session_find_halos.restype = ctypes.c_int
        lib.enzomodules_session_find_halos.argtypes = [
            ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int]
        lib.enzomodules_session_halo_count.restype = ctypes.c_int
        lib.enzomodules_session_halo_count.argtypes = []
        lib.enzomodules_session_subhalo_count.restype = ctypes.c_int
        lib.enzomodules_session_subhalo_count.argtypes = []
        lib.enzomodules_session_halo_property.restype = ctypes.c_double
        lib.enzomodules_session_halo_property.argtypes = [
            ctypes.c_int, ctypes.c_int]
        for fn in ("advance_time", "update_particles", "copy_baryon_to_old",
                   "clear_boundary_fluxes"):
            getattr(lib, "enzomodules_session_" + fn).argtypes = [ctypes.c_void_p, ctypes.c_int]
        lib._problem_set = True
    return lib


class GridView:
    """A read-only view of one initialized grid in the hierarchy."""

    def __init__(self, lib, handle, index):
        self._lib = lib
        self._h = handle
        self.index = index

    @property
    def rank(self) -> int:
        return self._lib.enzomodules_problem_grid_rank(self._h, self.index)

    @property
    def dims(self) -> List[int]:
        d = (ctypes.c_int * 3)()
        self._lib.enzomodules_problem_grid_dims(self._h, self.index, d)
        return [d[0], d[1], d[2]]

    @property
    def size(self) -> int:
        return self._lib.enzomodules_problem_grid_size(self._h, self.index)

    @property
    def field_types(self) -> List[int]:
        n = self._lib.enzomodules_problem_num_fields(self._h, self.index)
        t = (ctypes.c_int * n)()
        self._lib.enzomodules_problem_field_types(self._h, self.index, t)
        return [t[i] for i in range(n)]

    @property
    def field_names(self) -> List[str]:
        return [FIELD_NAMES.get(t, f"Field{t}") for t in self.field_types]

    def field(self, which) -> List[float]:
        """Field data (flat, includes ghost zones).  ``which`` is a field
        index, a FieldType int, or a name like ``"Density"``."""
        types = self.field_types
        if isinstance(which, str):
            fi = types.index(NAME_TO_FIELD[which])
        elif which in types and which not in range(len(types)):
            fi = types.index(which)
        else:
            fi = which
        out = (ctypes.c_double * self.size)()
        self._lib.enzomodules_problem_get_field(self._h, self.index, fi, out)
        return list(out)

    @property
    def num_particles(self) -> int:
        return self._lib.enzomodules_problem_num_particles(self._h, self.index)

    def particle_positions(self, dim: int) -> List[float]:
        n = self.num_particles
        out = (ctypes.c_double * n)()
        self._lib.enzomodules_problem_get_particle_pos(self._h, self.index, dim, out)
        return list(out)


class Problem:
    """Initialized Enzo problem (a grid hierarchy) built from a parameter file."""

    def __init__(self, paramfile: str, workdir: Optional[str] = None,
                 evolve: bool = False, stop_time: float = 0.0,
                 stop_cycle: int = 0):
        """Initialize a problem from ``paramfile``.  If ``evolve`` is True, also
        run Enzo's full time integration (EvolveHierarchy) to ``stop_time`` /
        ``stop_cycle`` (0 = use the parameter file's value); the resulting
        handle holds the *evolved* hierarchy."""
        if not available():
            raise RuntimeError(
                f"grid solver library not built ({bridge.grid_libpath()}); "
                "run EnzoModules/deps/build_grid.sh")
        self._lib = _lib()
        paramfile = os.path.abspath(paramfile)
        # InitializeNew / EvolveHierarchy write log + dump files to cwd; isolate.
        self._tmp = None
        cwd = os.getcwd()
        if workdir is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="enzomodules_ic_")
            workdir = self._tmp.name
        try:
            os.chdir(workdir)
            if evolve:
                with _suppress_fd_output():       # quiet the cycle log
                    self._h = self._lib.enzomodules_evolve_problem(
                        paramfile.encode(), float(stop_time), int(stop_cycle))
            else:
                self._h = self._lib.enzomodules_init_problem(paramfile.encode())
        finally:
            os.chdir(cwd)
        if not self._h:
            verb = "EvolveHierarchy" if evolve else "InitializeNew"
            raise RuntimeError(f"{verb} failed for {paramfile}")

    @property
    def problem_type(self) -> int:
        return self._lib.enzomodules_problem_problemtype(self._h)

    @property
    def num_grids(self) -> int:
        return self._lib.enzomodules_problem_num_grids(self._h)

    def grid(self, index: int = 0) -> GridView:
        return GridView(self._lib, self._h, index)

    def close(self):
        if getattr(self, "_h", None):
            self._lib.enzomodules_free_problem(self._h)
            self._h = None
        if self._tmp is not None:
            self._tmp.cleanup()
            self._tmp = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class Session:
    """A *live* Enzo hierarchy whose time loop you drive from Python.

    Where :class:`Problem` (``evolve=True``) hands the whole run to Enzo's
    monolithic ``EvolveHierarchy``, a ``Session`` exposes each orchestration
    step on the live hierarchy -- boundary conditions, the CFL timestep, the
    hydro/MHD solver, time advance and regridding -- so a host language can own
    the loop (e.g. a Python re-implementation of ``EvolveLevel``) while calling
    the certified legacy reference for each step.  Read the evolving state with
    the same :class:`GridView` accessors as :class:`Problem`.

        with Session("Toro-1-ShockTube.enzo") as s:
            s.run(level=0)             # default loop: == EvolveHierarchy
            rho = s.grid(0).field("Density")

    or step manually::

        while s.time < s.stop_time:
            s.set_boundary(0)
            dt = s.compute_dt(0); s.set_dt(0, dt)
            s.solve_hydro(0)
            s.advance_time(0)
    """

    def __init__(self, paramfile: str, workdir: Optional[str] = None,
                 stop_time: float = 0.0):
        if not available():
            raise RuntimeError(
                f"grid solver library not built ({bridge.grid_libpath()}); "
                "run EnzoModules/deps/build_grid.sh")
        self._lib = _lib()
        paramfile = os.path.abspath(paramfile)
        self._tmp = None
        cwd = os.getcwd()
        if workdir is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="enzomodules_sess_")
            workdir = self._tmp.name
        self._workdir = os.path.abspath(workdir)
        try:
            os.chdir(workdir)
            with _suppress_fd_output():
                self._h = self._lib.enzomodules_session_init(paramfile.encode())
        finally:
            os.chdir(cwd)
        if not self._h:
            raise RuntimeError(f"session_init (InitializeNew) failed for {paramfile}")
        if stop_time > 0.0:
            # Override StopTime by stepping nothing; reuse the C field directly
            # is not exposed, so callers pass stop_time to run() instead.
            self._stop_override = stop_time
        else:
            self._stop_override = None

    # --- live state -------------------------------------------------------
    @property
    def time(self) -> float:
        return self._lib.enzomodules_session_time(self._h)

    @property
    def stop_time(self) -> float:
        if self._stop_override is not None:
            return self._stop_override
        return self._lib.enzomodules_session_stop_time(self._h)

    @property
    def cycle(self) -> int:
        return self._lib.enzomodules_session_cycle(self._h)

    @property
    def num_grids(self) -> int:
        return self._lib.enzomodules_problem_num_grids(self._h)

    def grid(self, index: int = 0) -> GridView:
        return GridView(self._lib, self._h, index)

    # --- orchestration steps (one level at a time) ------------------------
    def set_boundary(self, level: int = 0) -> None:
        if self._lib.enzomodules_session_set_boundary(self._h, level):
            raise RuntimeError(f"SetBoundaryConditions failed (level {level})")

    def compute_dt(self, level: int = 0) -> float:
        return self._lib.enzomodules_session_compute_dt(self._h, level)

    def set_dt(self, level: int, dt: float) -> None:
        self._lib.enzomodules_session_set_dt(self._h, level, float(dt))

    def solve_hydro(self, level: int = 0) -> None:
        if self._lib.enzomodules_session_solve_hydro(self._h, level):
            raise RuntimeError(f"SolveHydroEquations failed (level {level})")

    def solve_cooling(self, level: int = 0) -> None:
        """Solve radiative cooling + non-equilibrium species rate equations
        (grid::MultiSpeciesHandler) -- the cooling/chemistry sub-step EvolveLevel
        runs after the hydro solve.  Dispatches to Grackle, the coupled
        rate-and-cool solver, or SolveRateEquations + SolveRadiativeCooling per
        the run's settings.  No-op unless MultiSpecies or RadiativeCooling is on.
        Call after solve_hydro(level)."""
        if self._lib.enzomodules_session_solve_cooling(self._h, level):
            raise RuntimeError(f"cooling/chemistry solve failed (level {level})")

    def star_particles(self, level: int = 0) -> None:
        """Run the star-particle lifecycle on `level` (the same sequence
        EvolveLevel drives): StarParticleInitialize (gather stars, set feedback
        flags) -> per-grid StarParticleHandler (star formation + feedback) ->
        StarParticleFinalize (ActivateNewStar -- flip newly-formed/aged particles
        to active, radiating stars, and sync back to the grids).  No-op unless
        StarParticleCreation or StarParticleFeedback is on.  Call after
        solve_hydro/solve_cooling; a later evolve_photons(stars=True) then
        radiates from the now-active stars.

        Note: a star only enters its radiating main-sequence window for its
        (often short, e.g. Pop III) lifetime, so resolving it requires timesteps
        finer than that lifetime."""
        if self._lib.enzomodules_session_star_particles(self._h, level):
            raise RuntimeError(f"star-particle lifecycle failed (level {level})")

    def advance_time(self, level: int = 0) -> None:
        self._lib.enzomodules_session_advance_time(self._h, level)

    def rebuild(self, level: int = 0) -> None:
        if self._lib.enzomodules_session_rebuild(self._h, level):
            raise RuntimeError(f"RebuildHierarchy failed (level {level})")

    def gravity(self, level: int = 0) -> None:
        """Self-gravity chain: deposit mass, solve Poisson (PrepareDensityField),
        then per-grid accelerations + potential.  No-op unless SelfGravity is on.
        Call after set_boundary(level)."""
        if self._lib.enzomodules_session_gravity(self._h, level):
            raise RuntimeError(f"gravity chain failed (level {level})")

    def update_particles(self, level: int = 0) -> None:
        """Drift particle positions by their grid's dtFixed."""
        self._lib.enzomodules_session_update_particles(self._h, level)

    def evolve_photons(self, level: int = 0, stars: bool = False) -> None:
        """Emit and transport radiative-transfer photon packages, depositing
        photo-ionization / heating rates and coupling to the chemistry.  Sets up
        the SubgridMarker grid-ownership map, sizes the photon timestep, and
        sub-cycles the photon time up to the grid time so the sources radiate.

        ``stars=False`` (default) uses the radiation sources read from the
        parameter file (e.g. PhotonTest).  ``stars=True`` instead gathers the
        grid's star particles (StarParticleInitialize) and converts the active
        ones into radiation sources (StarParticleRadTransfer), so star particles
        radiate -- note a star only emits once Enzo's formation lifecycle has
        marked it an active radiation source (type > 0, alive, non-formation
        feedback flag).

        No-op unless RadiativeTransfer is on (library built with -DTRANSFER).
        Call after set_dt(level, ...)."""
        if self._lib.enzomodules_session_evolve_photons_ex(
                self._h, level, 1 if stars else 0):
            raise RuntimeError(f"EvolvePhotons failed (level {level})")

    def find_halos(self, subfind: bool = False, linking_length: float = 0.0,
                   min_size: int = 0) -> list:
        """Run Enzo's inline friends-of-friends halo finder (and SUBFIND if
        ``subfind=True``) on the session's particles, returning a list of halo
        dicts sorted largest-first.  Each dict has: ``n`` (particle count),
        ``mass`` / ``virial_mass`` [Msun], ``virial_radius`` [kpc], ``center``
        (x,y,z centre of mass), ``velocity`` (mean, km/s), ``vrms`` [km/s],
        ``spin``, ``angular_momentum`` (x,y,z).  The total number of SUBFIND
        subgroups is available afterwards via ``subhalo_count``.

        ``linking_length`` (FOF b, in mean-interparticle-separation units) and
        ``min_size`` (minimum particles per halo; clamped to >= 32, SUBFIND's
        neighbour requirement) override HaloFinderLinkingLength /
        HaloFinderMinimumSize when > 0.

        The finder moves the particles off the grids and restores them
        (FOF_Initialize/FOF_Finalize), so this is repeatable and leaves the
        hierarchy's particle distribution intact (total count conserved).
        Returns [] if the build has no particles."""
        n = self._lib.enzomodules_session_find_halos(
            self._h, 1 if subfind else 0, float(linking_length), int(min_size))
        if n < 0:
            raise RuntimeError("inline halo finder failed")
        P = self._lib.enzomodules_session_halo_property
        halos = []
        for i in range(n):
            halos.append({
                "n":                int(P(i, 0)),
                "mass":             P(i, 1),
                "center":           (P(i, 2), P(i, 3), P(i, 4)),
                "velocity":         (P(i, 5), P(i, 6), P(i, 7)),
                "vrms":             P(i, 8),
                "spin":             P(i, 9),
                "virial_mass":      P(i, 10),
                "virial_radius":    P(i, 11),
                "angular_momentum": (P(i, 12), P(i, 13), P(i, 14)),
            })
        return halos

    def subhalo_count(self) -> int:
        """Total number of SUBFIND subgroups from the most recent
        ``find_halos(subfind=True)`` (0 if SUBFIND was not run)."""
        return self._lib.enzomodules_session_subhalo_count()

    # --- AMR conservation machinery (for a multi-level EvolveLevel) --------
    def num_grids_on_level(self, level: int) -> int:
        """Number of grids on `level` (0 if the level is empty)."""
        return self._lib.enzomodules_session_num_grids_on_level(self._h, level)

    def create_fluxes(self, level: int = 0) -> None:
        """Allocate the per-level boundary-flux storage that solve_hydro fills
        and update_from_finer reads.  Call before solve_hydro when you need
        conservative AMR flux correction."""
        if self._lib.enzomodules_session_create_fluxes(self._h, level):
            raise RuntimeError(f"CreateFluxes failed (level {level})")

    def clear_boundary_fluxes(self, level: int = 0) -> None:
        """Zero the boundary-flux accumulators on all grids of `level`."""
        self._lib.enzomodules_session_clear_boundary_fluxes(self._h, level)

    def copy_baryon_to_old(self, level: int = 0) -> None:
        """Save the current baryon fields as the 'old' fields (time-centering)."""
        self._lib.enzomodules_session_copy_baryon_to_old(self._h, level)

    def update_from_finer(self, level: int = 0) -> None:
        """Conservative fine->coarse coupling: project the finer level into this
        one and correct boundary zones for the flux difference.  Requires
        create_fluxes + solve_hydro on this level and an evolved level+1."""
        if self._lib.enzomodules_session_update_from_finer(self._h, level):
            raise RuntimeError(f"UpdateFromFinerGrids failed (level {level})")

    def finalize_fluxes(self, level: int = 0) -> None:
        """Release the per-level flux storage allocated by create_fluxes."""
        if self._lib.enzomodules_session_finalize_fluxes(self._h, level):
            raise RuntimeError(f"FinalizeFluxes failed (level {level})")

    def evolve_level(self, level: int = 0, dt_above: float = 0.0,
                     gravity: bool = False, regrid: bool = True,
                     radiation: bool = False, star_sources: bool = False,
                     cooling: bool = False, star_formation: bool = False) -> int:
        """A Python re-implementation of Enzo's recursive EvolveLevel, built
        entirely from the certified session steps.  Mirrors the legacy control
        flow: clear the boundary fluxes once on entry, then sub-cycle this level
        until it catches up to its parent's timestep (``dt_above``) -- each
        sub-cycle solving the grids, recursing into level+1, then applying
        conservative flux correction + projection (``update_from_finer``) on the
        way back up, and regridding the finer levels between sub-cycles.
        ``dt_above == 0`` means a single (top-grid) step.  Returns the number of
        sub-cycles taken.

        This is the proof that EvolveLevel *can* be written in Python on top of
        these bridges -- the full AMR time integrator, step by certified step."""
        self.clear_boundary_fluxes(level)
        n = 0
        done = 0.0
        while True:
            self.set_boundary(level)                # interpolate from parent
            dt = self.compute_dt(level)
            if dt_above > 0.0:
                dt = min(dt, dt_above - done)        # don't overshoot the parent
            self.set_dt(level, dt)

            if radiation:                           # emit + transport photons
                self.evolve_photons(level, stars=star_sources)   # -> chemistry

            self.create_fluxes(level)               # allocate flux storage
            if gravity:
                self.gravity(level)
            self.copy_baryon_to_old(level)
            self.solve_hydro(level)                 # fills the boundary fluxes
            if cooling:
                self.solve_cooling(level)           # radiative cooling + chemistry
            if star_formation:
                self.star_particles(level)          # form / feedback / activate
            self.update_particles(level)
            self.advance_time(level)

            last = (dt_above <= 0.0) or (done + dt >= dt_above * (1 - 1e-6))

            if self.num_grids_on_level(level + 1) > 0:
                self.set_boundary(level)            # refresh before projection
                self.evolve_level(level + 1, dt_above=dt, gravity=gravity,
                                  regrid=regrid, radiation=radiation,
                                  star_sources=star_sources, cooling=cooling,
                                  star_formation=star_formation)
                self.update_from_finer(level)       # project + flux-correct
            self.finalize_fluxes(level)

            n += 1
            done += dt
            if last or dt <= 0:
                break
            if regrid:
                self.rebuild(level)                 # regrid finer levels
        return n

    def run_amr(self, gravity: bool = False, regrid: bool = True,
                radiation: bool = False, star_sources: bool = False,
                cooling: bool = False, star_formation: bool = False,
                max_cycles: int = 100000) -> int:
        """Top-level driver mirroring EvolveHierarchy: an initial regrid, then
        repeatedly evolve the whole hierarchy one root step (``evolve_level(0)``)
        and regrid, until StopTime.  Pass ``gravity=True`` to include the
        self-gravity chain, ``radiation=True`` to emit + transport photons each
        step, ``star_sources=True`` to radiate from star particles,
        ``cooling=True`` to solve radiative cooling + chemistry, and
        ``star_formation=True`` to run the star formation/feedback/activation
        lifecycle each step.  This is a complete (optionally
        radiation-hydrodynamic + gravitating + cooling + star-forming) AMR run
        driven from Python on top of the certified legacy steps."""
        n = 0
        with _suppress_fd_output():
            if regrid:
                self.rebuild(0)
            while self.time < self.stop_time and n < max_cycles:
                self.evolve_level(0, gravity=gravity, regrid=regrid,
                                  radiation=radiation, star_sources=star_sources,
                                  cooling=cooling, star_formation=star_formation)
                if regrid:
                    self.rebuild(0)
                n += 1
        return n

    # --- the default top-level loop (reproduces EvolveHierarchy on level 0)-
    def step(self, level: int = 0) -> float:
        """One root-level step: boundary -> CFL dt -> solve -> advance.
        Returns the dt taken."""
        self.set_boundary(level)
        dt = self.compute_dt(level)
        self.set_dt(level, dt)
        self.solve_hydro(level)
        self.advance_time(level)
        return dt

    def run(self, level: int = 0, max_cycles: int = 100000) -> int:
        """Drive the loop to StopTime.  Returns the number of cycles taken."""
        n = 0
        with _suppress_fd_output():
            while self.time < self.stop_time and n < max_cycles:
                dt = self.step(level)
                if dt <= 0:
                    break
                n += 1
        return n

    # --- data output / restart -------------------------------------------
    def write_output(self, number: int = 0, checkpoint: bool = False) -> str:
        """Write the full simulation state to disk (Group_WriteAllData): the
        same HDF5 data dump Enzo writes -- grid data, hierarchy, external
        boundary and a parameter file -- into this session's working directory
        (named from the run's DataDumpName + ``number``).  ``checkpoint=True``
        flags it as a checkpoint.  Returns the absolute path of the dump's
        hierarchy/parameter file, ready to pass to ``Session.from_output``."""
        cwd = os.getcwd()
        os.chdir(self._workdir)
        try:
            with _suppress_fd_output():
                rc = self._lib.enzomodules_session_write_output(
                    self._h, int(number), 1 if checkpoint else 0)
            if rc:
                raise RuntimeError("write_output (Group_WriteAllData) failed")
            # Locate the dump by its .hierarchy file (robust to Enzo's naming).
            hier = [f for f in os.listdir(".") if f.endswith(".hierarchy")]
            if not hier:
                raise RuntimeError("write_output produced no hierarchy file")
            hier.sort(key=lambda f: os.path.getmtime(f))   # newest = this dump
            name = hier[-1][:-len(".hierarchy")]
            return os.path.join(self._workdir, name)
        finally:
            os.chdir(cwd)

    @classmethod
    def from_output(cls, path: str) -> "Session":
        """Reload a session from a dump written by ``write_output`` (or by Enzo
        itself): Group_ReadAllData restores the parameter file, external
        boundary, hierarchy and grid data.  ``path`` is the dump's
        hierarchy/parameter file (the value ``write_output`` returned, with or
        without the ``.hierarchy`` suffix)."""
        if not available():
            raise RuntimeError("grid solver library not built")
        path = os.path.abspath(path)
        if path.endswith(".hierarchy"):
            path = path[:-len(".hierarchy")]
        workdir = os.path.dirname(path)
        name = os.path.basename(path)
        self = cls.__new__(cls)
        self._lib = _lib()
        self._tmp = None
        self._workdir = workdir
        self._stop_override = None
        cwd = os.getcwd()
        os.chdir(workdir)
        try:
            with _suppress_fd_output():
                self._h = self._lib.enzomodules_session_from_output(name.encode())
        finally:
            os.chdir(cwd)
        if not self._h:
            raise RuntimeError(f"from_output (Group_ReadAllData) failed for {path}")
        return self

    def close(self):
        if getattr(self, "_h", None):
            self._lib.enzomodules_free_problem(self._h)
            self._h = None
        if self._tmp is not None:
            self._tmp.cleanup()
            self._tmp = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
