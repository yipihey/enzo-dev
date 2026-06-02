"""A 1D finite-volume hydro driver built on the wrapped hydro_rk Riemann solver.

Every timestep computes interface fluxes with
``enzomodules.bridge.hydro_rk_line`` (the legacy Enzo HLL/HLLC/LLF + PLM line
solver) and does a conservative update.  This is the hydro_rk analogue of
``ppm_sod`` and the reference for orchestrating the wrapped C++ solver.

Conserved state U = [rho, rho*vx, rho*vy, rho*vz, rho*Etot] (Enzo flux index
order iD,iS1,iS2,iS3,iEtot); primitives prim = [rho, eint, vx, vy, vz].
Time integration is forward-Euler MUSCL (PLM in space); small CFL keeps it
accurate enough to validate against the exact Riemann solution.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from .. import bridge

NEQ = 5


def _primitives(U: List[List[float]], gamma: float) -> List[List[float]]:
    n = len(U[0])
    rho = U[0]
    vx = [U[1][i] / rho[i] for i in range(n)]
    vy = [U[2][i] / rho[i] for i in range(n)]
    vz = [U[3][i] / rho[i] for i in range(n)]
    etot = [U[4][i] / rho[i] for i in range(n)]
    eint = [etot[i] - 0.5 * (vx[i] ** 2 + vy[i] ** 2 + vz[i] ** 2) for i in range(n)]
    return [list(rho), eint, vx, vy, vz]


def _fill_ghosts(U: List[List[float]], nghost: int) -> None:
    n = len(U[0])
    lo, hi = nghost, n - 1 - nghost
    for f in range(NEQ):
        for k in range(nghost):
            U[f][k] = U[f][lo]
            U[f][n - 1 - k] = U[f][hi]


def _max_speed(prim, gamma, nghost):
    rho, eint, vx = prim[0], prim[1], prim[2]
    n = len(rho)
    cmax = 0.0
    for c in range(nghost, n - nghost):
        p = (gamma - 1.0) * rho[c] * eint[c]
        cs = math.sqrt(gamma * max(p, 1e-20) / rho[c])
        cmax = max(cmax, abs(vx[c]) + cs)
    return cmax


def sod_state(nx, nghost, gamma, left, right, disc, length):
    dx = length / nx
    ncells = nx + 2 * nghost
    U = [[0.0] * ncells for _ in range(NEQ)]
    for c in range(ncells):
        x = (c - nghost + 0.5) * dx
        rho, u, p = left if x < disc else right
        etot = p / ((gamma - 1.0) * rho) + 0.5 * u * u
        U[0][c] = rho
        U[1][c] = rho * u
        U[2][c] = 0.0
        U[3][c] = 0.0
        U[4][c] = rho * etot
    return U, dx


def run(solver="hll", t_final=0.2, nx=200, cfl=0.3, gamma=1.4, nghost=3,
        theta_limiter=1.5, length=1.0, disc=0.5,
        left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1),
        max_steps=100000) -> Tuple[List[float], List[float], List[float], List[float], int]:
    """Evolve a Sod tube to ``t_final``.  Returns (x, rho, u, p, nsteps)."""
    code = {"hll": bridge.HLL, "llf": bridge.LLF, "hllc": bridge.HLLC}[solver]
    U, dx = sod_state(nx, nghost, gamma, left, right, disc, length)
    t = 0.0
    n = 0
    while t < t_final and n < max_steps:
        _fill_ghosts(U, nghost)
        prim = _primitives(U, gamma)
        dt = cfl * dx / _max_speed(prim, gamma, nghost)
        if t + dt > t_final:
            dt = t_final - t
        flux = bridge.hydro_rk_line(code, prim, gamma=gamma,
                                    theta_limiter=theta_limiter, nghost=nghost)
        r = dt / dx
        for f in range(NEQ):
            for j in range(nx):
                c = nghost + j
                U[f][c] -= r * (flux[f][j + 1] - flux[f][j])
        t += dt
        n += 1

    x = [(j + 0.5) * dx for j in range(nx)]
    rho = [U[0][nghost + j] for j in range(nx)]
    u = [U[1][nghost + j] / rho[j] for j in range(nx)]
    p = [(gamma - 1.0) * rho[j] * (U[4][nghost + j] / rho[j] - 0.5 * u[j] ** 2)
         for j in range(nx)]
    return x, rho, u, p, n
