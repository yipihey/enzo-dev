"""Ergonomic wrappers for the legacy Enzo hydro kernels.

One module per kernel group; this is the first (``twoshock.F``).
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from . import bridge
from .fixtures import Fixture


def twoshock(dls: Sequence[float], drs: Sequence[float],
             pls: Sequence[float], prs: Sequence[float],
             uls: Sequence[float], urs: Sequence[float],
             *, gamma: float,
             dt: float = 0.0, pmin: float = 1e-20,
             i1: int = 1, i2: Optional[int] = None,
             j1: int = 1, j2: int = 1,
             ipresfree: int = 0, gravity: int = 0,
             idual: int = 0, eta1: float = 0.0,
             grslice: Optional[Sequence[float]] = None
             ) -> Tuple[List[float], List[float]]:
    """Two-shock approximate Riemann solver (legacy ``twoshock.F``).

    Inputs are the reconstructed left/right interface states; returns the
    resolved interface pressure ``pbar`` and normal velocity ``ubar`` as
    lists of length ``len(dls)``.  Single-slice convenience: pass length-N
    sequences and the full range is solved by default.
    """
    idim = len(dls)
    jdim = 1
    if i2 is None:
        i2 = idim
    if grslice is None:
        grslice = [0.0] * idim
    return bridge.twoshock_raw(
        dls, drs, pls, prs, uls, urs,
        idim, jdim, i1, i2, j1, j2,
        dt, gamma, pmin, ipresfree,
        gravity, grslice, idual, eta1,
    )


def ppm_sweep_1d(dslice, eslice, uslice, vslice, wslice, pslice,
                 *, i1: int, i2: int, dx: float, dt: float, gamma: float,
                 want_fluxes: bool = False):
    """One PPM hydro update of a 1D slice (legacy ``Grid_xEulerSweep`` core).

    Composes ``inteuler -> twoshock -> flux_twoshock -> euler`` for the
    simplest faithful configuration (no gravity / dual-energy / colour /
    diffusion; PPM + two-shock).  ``dslice`` (density), ``eslice`` (total
    specific energy), ``uslice``/``vslice``/``wslice`` (velocities) and
    ``pslice`` (pressure, precomputed by the caller) are length-``idim``
    sequences with >= 3 ghost cells each side; ``i1``/``i2`` are the 1-based
    inclusive active range.

    Returns ``(d, e, u, v, w)`` updated lists, plus ``(df, ef, uf)`` if
    ``want_fluxes`` else ``None``.  See ``enzomodules.examples.ppm_sod`` for a
    full evolution driver built on this.
    """
    return bridge.ppm_sweep_1d(dslice, eslice, uslice, vslice, wslice, pslice,
                               i1, i2, dx, dt, gamma, want_fluxes=want_fluxes)


def ppm_sweep_from_fixture(fx: Fixture):
    """Run one PPM sweep using a fixture's stored input slice, reproducing how
    its reference outputs were captured.  Returns ``(d, e, u, v, w)``."""
    state, _ = bridge.ppm_sweep_1d(
        fx["dslice_in"], fx["eslice_in"], fx["uslice_in"],
        fx["vslice_in"], fx["wslice_in"], fx["pslice_in"],
        int(fx["i1"]), int(fx["i2"]),
        float(fx["dx"]), float(fx["dt"]), float(fx["gamma"]),
        want_fluxes=False)
    return state


def twoshock_from_fixture(fx: Fixture) -> Tuple[List[float], List[float]]:
    """Run the kernel using a fixture's stored inputs, reproducing exactly
    how its reference outputs were captured.  Returns the freshly computed
    ``(pbar, ubar)`` for comparison against ``fx['pbar']`` / ``fx['ubar']``."""
    return bridge.twoshock_raw(
        fx["dls"], fx["drs"], fx["pls"], fx["prs"], fx["uls"], fx["urs"],
        int(fx["idim"]), int(fx["jdim"]),
        int(fx["i1"]), int(fx["i2"]), int(fx["j1"]), int(fx["j2"]),
        float(fx["dt"]), float(fx["gamma"]), float(fx["pmin"]),
        int(fx["ipresfree"]),
        int(fx["gravity"]), fx["grslice"], int(fx["idual"]),
        float(fx["eta1"]),
    )
