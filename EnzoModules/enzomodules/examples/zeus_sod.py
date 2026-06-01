"""1D hydro driver built on the wrapped ZEUS solver (grid::ZeusSolver).

Demonstrates the grid-fixture path: every step calls
``enzomodules.bridge.zeus_sweep_1d``, which builds a minimal Enzo grid, runs
the operator-split finite-difference ZEUS update, and reads the fields back.

ZEUS is internal-energy based, so the energy field holds specific internal
energy e = p / ((gamma-1) rho) (no kinetic part).  State is cell-centred
[rho, e, vx]; forward time-stepping with zero-gradient (outflow) boundaries.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from .. import bridge

NGHOST = 3


def sod_state(nx, gamma, length, disc, left, right):
    dx = length / nx
    idim = nx + 2 * NGHOST
    d = [0.0] * idim
    e = [0.0] * idim
    u = [0.0] * idim
    for c in range(idim):
        x = (c - NGHOST + 0.5) * dx
        rho, vel, p = left if x < disc else right
        d[c] = rho
        u[c] = vel
        e[c] = p / ((gamma - 1.0) * rho)      # specific internal energy
    return d, e, u, dx


def _fill_ghosts(d, e, u):
    n = len(d)
    lo, hi = NGHOST, n - 1 - NGHOST
    for a in (d, e, u):
        for k in range(NGHOST):
            a[k] = a[lo]
            a[n - 1 - k] = a[hi]


def _max_speed(d, e, u, gamma):
    n = len(d)
    cmax = 0.0
    for c in range(NGHOST, n - NGHOST):
        p = (gamma - 1.0) * d[c] * e[c]
        cs = math.sqrt(gamma * max(p, 1e-20) / d[c])
        cmax = max(cmax, abs(u[c]) + cs)
    return cmax


def run(t_final=0.2, nx=200, cfl=0.4, gamma=1.4, length=1.0, disc=0.5,
        left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1),
        max_steps=100000) -> Tuple[List[float], List[float], List[float], List[float], int]:
    """Evolve a Sod tube with ZEUS to ``t_final``.  Returns (x, rho, u, p, n)."""
    d, e, u, dx = sod_state(nx, gamma, length, disc, left, right)
    t = 0.0
    n = 0
    while t < t_final and n < max_steps:
        _fill_ghosts(d, e, u)
        dt = cfl * dx / _max_speed(d, e, u, gamma)
        if t + dt > t_final:
            dt = t_final - t
        d, e, u = bridge.zeus_sweep_1d(d, e, u, NGHOST, dx, dt, gamma)
        t += dt
        n += 1

    x = [(j + 0.5) * dx for j in range(nx)]
    rho = [d[NGHOST + j] for j in range(nx)]
    vel = [u[NGHOST + j] for j in range(nx)]
    p = [(gamma - 1.0) * d[NGHOST + j] * e[NGHOST + j] for j in range(nx)]
    return x, rho, vel, p, n
