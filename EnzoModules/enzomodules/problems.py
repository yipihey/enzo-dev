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


def available() -> bool:
    return bridge.grid_available()


def _lib():
    lib = bridge._load_gridlib()
    if not getattr(lib, "_problem_set", False):
        d = ctypes.POINTER(ctypes.c_double)
        ip = ctypes.POINTER(ctypes.c_int)
        lib.enzomodules_init_problem.restype = ctypes.c_void_p
        lib.enzomodules_init_problem.argtypes = [ctypes.c_char_p]
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

    def __init__(self, paramfile: str, workdir: Optional[str] = None):
        if not available():
            raise RuntimeError(
                f"grid solver library not built ({bridge.grid_libpath()}); "
                "run EnzoModules/deps/build_grid.sh")
        self._lib = _lib()
        paramfile = os.path.abspath(paramfile)
        # InitializeNew writes OutputLog etc. to cwd; isolate it.
        self._tmp = None
        cwd = os.getcwd()
        if workdir is None:
            self._tmp = tempfile.TemporaryDirectory(prefix="enzomodules_ic_")
            workdir = self._tmp.name
        try:
            os.chdir(workdir)
            self._h = self._lib.enzomodules_init_problem(paramfile.encode())
        finally:
            os.chdir(cwd)
        if not self._h:
            raise RuntimeError(f"InitializeNew failed for {paramfile}")

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
