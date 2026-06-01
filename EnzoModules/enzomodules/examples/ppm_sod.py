"""A complete 1D PPM hydro driver built on a single certified kernel.

This is the reference for how a downstream rewrite orchestrates a wrapped
solver: everything below (initial conditions, boundary fill, pressure,
CFL-limited timestep, time loop) is ordinary host-language code, and the only
physics call is ``enzomodules.hydro.ppm_sweep_1d`` -- the legacy Enzo PPM
sweep, certified against golden fixtures.

The structure mirrors Enzo's own ``Grid_SolvePPM_DE`` -> ``xEulerSweep``:
compute pressure once per step, refresh boundaries, advance one sweep.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple

from .. import hydro

NGHOST = 3  # PPM stencil needs 3 ghost cells per side


@dataclass
class Grid1D:
    """A 1D hydro state on a uniform grid with ghost zones."""
    nx: int
    dx: float
    gamma: float
    d: List[float]
    e: List[float]   # total specific energy
    u: List[float]
    v: List[float]
    w: List[float]
    x0: float = 0.0  # left edge of the active domain

    @property
    def idim(self) -> int:
        return self.nx + 2 * NGHOST

    @property
    def i1(self) -> int:           # 1-based first active cell
        return NGHOST + 1

    @property
    def i2(self) -> int:           # 1-based last active cell
        return NGHOST + self.nx

    def x_centers(self) -> List[float]:
        """Cell-center positions of the active cells."""
        return [self.x0 + (i + 0.5) * self.dx for i in range(self.nx)]

    def pressure(self) -> List[float]:
        g = self.gamma
        return [(g - 1.0) * self.d[c] * (self.e[c] - 0.5 * self.u[c] ** 2)
                for c in range(self.idim)]

    def active(self, arr: List[float]) -> List[float]:
        return arr[NGHOST:self.idim - NGHOST]


def sod_grid(nx: int = 200, length: float = 1.0, gamma: float = 1.4,
             discontinuity: float = 0.5,
             left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1)) -> Grid1D:
    """Construct a Sod shock-tube initial state."""
    dx = length / nx
    idim = nx + 2 * NGHOST
    dL, uL, pL = left
    dR, uR, pR = right
    d = [0.0] * idim
    u = [0.0] * idim
    e = [0.0] * idim
    for c in range(idim):
        xc = (c - NGHOST + 0.5) * dx
        if xc < discontinuity:
            dc, uc, pc = dL, uL, pL
        else:
            dc, uc, pc = dR, uR, pR
        d[c], u[c] = dc, uc
        e[c] = pc / ((gamma - 1.0) * dc) + 0.5 * uc * uc
    return Grid1D(nx=nx, dx=dx, gamma=gamma, d=d, e=e, u=u,
                  v=[0.0] * idim, w=[0.0] * idim, x0=0.0)


def _fill_ghosts(g: Grid1D) -> None:
    """Zero-gradient (outflow) boundaries."""
    lo, hi = NGHOST, g.idim - 1 - NGHOST
    for k in range(NGHOST):
        for a in (g.d, g.e, g.u, g.v, g.w):
            a[k] = a[lo]
            a[g.idim - 1 - k] = a[hi]


def max_wave_speed(g: Grid1D, p: List[float]) -> float:
    cmax = 0.0
    for c in range(NGHOST, g.idim - NGHOST):
        cs = math.sqrt(g.gamma * max(p[c], 1e-20) / g.d[c])
        cmax = max(cmax, abs(g.u[c]) + cs)
    return cmax


def step(g: Grid1D, dt: float) -> None:
    """Advance the grid by one PPM sweep of size ``dt`` (in place)."""
    p = g.pressure()
    (g.d, g.e, g.u, g.v, g.w), _ = hydro.ppm_sweep_1d(
        g.d, g.e, g.u, g.v, g.w, p,
        i1=g.i1, i2=g.i2, dx=g.dx, dt=dt, gamma=g.gamma)


def run(t_final: float, nx: int = 200, cfl: float = 0.4,
        gamma: float = 1.4, max_steps: int = 100000,
        **ic) -> Tuple[Grid1D, int]:
    """Evolve a Sod shock tube to ``t_final``; returns (grid, nsteps)."""
    g = sod_grid(nx=nx, gamma=gamma, **ic)
    t = 0.0
    n = 0
    while t < t_final and n < max_steps:
        _fill_ghosts(g)
        p = g.pressure()
        dt = cfl * g.dx / max_wave_speed(g, p)
        if t + dt > t_final:
            dt = t_final - t
        step(g, dt)
        t += dt
        n += 1
    return g, n
