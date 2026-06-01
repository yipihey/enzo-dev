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
        for fn in ("time", "stop_time", "compute_dt"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = ctypes.c_double
        lib.enzomodules_session_time.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_stop_time.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_cycle.restype = ctypes.c_int
        lib.enzomodules_session_cycle.argtypes = [ctypes.c_void_p]
        lib.enzomodules_session_compute_dt.argtypes = [ctypes.c_void_p, ctypes.c_int]
        for fn in ("set_boundary", "solve_hydro", "rebuild", "gravity",
                   "evolve_photons"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = ctypes.c_int
            f.argtypes = [ctypes.c_void_p, ctypes.c_int]
        for fn in ("set_dt", "advance_time", "update_particles"):
            f = getattr(lib, "enzomodules_session_" + fn)
            f.restype = None
        lib.enzomodules_session_set_dt.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double]
        lib.enzomodules_session_advance_time.argtypes = [ctypes.c_void_p, ctypes.c_int]
        lib.enzomodules_session_update_particles.argtypes = [ctypes.c_void_p, ctypes.c_int]
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

    def evolve_photons(self, level: int = 0) -> None:
        """Trace radiative-transfer photon packages.  No-op unless
        RadiativeTransfer is on (library built with -DTRANSFER)."""
        if self._lib.enzomodules_session_evolve_photons(self._h, level):
            raise RuntimeError(f"EvolvePhotons failed (level {level})")

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
