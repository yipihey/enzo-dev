"""1D MHD driver built on the wrapped Dedner-cleaning hydro_rk solver.

Evolves the Brio & Wu (1988) magnetized shock tube using
``enzomodules.bridge.mhd_rk_line`` (the legacy Enzo HLLD/HLL/LLF + PLM MHD
line solver with GLM divergence cleaning) for the interface fluxes, then a
conservative update.  The transverse-field flux carries the Dedner cleaning
(C_h), so an initially divergence-free Bx stays constant.

Conserved U = [rho, rho*vx, rho*vy, rho*vz, Etot, Bx, By, Bz, Phi] with
Etot = rho*(eint + 0.5 v^2) + 0.5 B^2; prim = [rho, eint, vx, vy, vz, Bx, By,
Bz, Phi].  Forward-Euler MUSCL in time (hyperbolic GLM only).
"""
from __future__ import annotations

import math
from typing import List, Tuple

from .. import bridge

NEQ = 9
BRIO_WU_LEFT = (1.0, 1.0, 1.0)     # rho, p, By
BRIO_WU_RIGHT = (0.125, 0.1, -1.0)
BRIO_WU_BX = 0.75
BRIO_WU_GAMMA = 2.0


def _primitives(U, gamma):
    n = len(U[0])
    P = [[0.0] * n for _ in range(NEQ)]
    for c in range(n):
        rho = U[0][c]
        vx, vy, vz = U[1][c] / rho, U[2][c] / rho, U[3][c] / rho
        Bx, By, Bz, Phi = U[5][c], U[6][c], U[7][c], U[8][c]
        B2 = Bx * Bx + By * By + Bz * Bz
        v2 = vx * vx + vy * vy + vz * vz
        eint = (U[4][c] - 0.5 * B2) / rho - 0.5 * v2
        P[0][c], P[1][c] = rho, eint
        P[2][c], P[3][c], P[4][c] = vx, vy, vz
        P[5][c], P[6][c], P[7][c], P[8][c] = Bx, By, Bz, Phi
    return P


def _fast_speed(rho, p, B2, Bx, gamma):
    gp = gamma * max(p, 1e-20)
    disc = max((gp + B2) ** 2 - 4.0 * gp * Bx * Bx, 0.0)
    return math.sqrt((gp + B2 + math.sqrt(disc)) / (2.0 * rho))


def _ch(P, gamma, nghost):
    n = len(P[0])
    ch = 0.0
    for c in range(nghost, n - nghost):
        rho = P[0][c]
        p = (gamma - 1.0) * rho * P[1][c]
        B2 = P[5][c] ** 2 + P[6][c] ** 2 + P[7][c] ** 2
        ch = max(ch, abs(P[2][c]) + _fast_speed(rho, p, B2, P[5][c], gamma))
    return ch


def brio_wu_state(nx, nghost, gamma, length, disc):
    dx = length / nx
    ncells = nx + 2 * nghost
    U = [[0.0] * ncells for _ in range(NEQ)]
    for c in range(ncells):
        x = (c - nghost + 0.5) * dx
        rho, p, By = BRIO_WU_LEFT if x < disc else BRIO_WU_RIGHT
        B2 = BRIO_WU_BX ** 2 + By ** 2
        U[0][c] = rho
        U[4][c] = p / (gamma - 1.0) + 0.5 * B2     # rho*eint + 0.5 B^2
        U[5][c] = BRIO_WU_BX
        U[6][c] = By
    return U, dx


def run(solver="hlld", t_final=0.1, nx=400, cfl=0.3, gamma=BRIO_WU_GAMMA,
        nghost=3, theta_limiter=1.5, length=1.0, disc=0.5, max_steps=100000):
    """Evolve Brio-Wu to ``t_final``.  Returns (x, rho, vx, By, Bx, nsteps)."""
    code = {"hlld": bridge.HLLD, "hll": bridge.HLL, "llf": bridge.LLF}[solver]
    U, dx = brio_wu_state(nx, nghost, gamma, length, disc)
    t = 0.0
    n = 0
    while t < t_final and n < max_steps:
        for f in range(NEQ):
            for k in range(nghost):
                U[f][k] = U[f][nghost]
                U[f][len(U[f]) - 1 - k] = U[f][len(U[f]) - 1 - nghost]
        P = _primitives(U, gamma)
        ch = _ch(P, gamma, nghost)
        dt = min(cfl * dx / ch, t_final - t)
        flux = bridge.mhd_rk_line(code, P, ch, gamma=gamma,
                                  theta_limiter=theta_limiter, nghost=nghost)
        r = dt / dx
        for f in range(NEQ):
            for j in range(nx):
                U[f][nghost + j] -= r * (flux[f][j + 1] - flux[f][j])
        t += dt
        n += 1

    x = [(j + 0.5) * dx for j in range(nx)]
    rho = [U[0][nghost + j] for j in range(nx)]
    vx = [U[1][nghost + j] / rho[j] for j in range(nx)]
    By = [U[6][nghost + j] for j in range(nx)]
    Bx = [U[5][nghost + j] for j in range(nx)]
    return x, rho, vx, By, Bx, n
