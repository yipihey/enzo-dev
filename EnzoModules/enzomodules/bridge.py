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
