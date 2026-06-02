"""2D hydro driver built on the wrapped full PPM step (grid::SolveHydroEquations).

Every step calls ``enzomodules.bridge.ppm_hydro_step``, which builds a minimal
2D Enzo grid, takes one PPM_DirectEuler step -- running BOTH the x and y
directional Euler sweeps (Grid_xEulerSweep / Grid_yEulerSweep) -- and reads the
fields back.

The grid is row-major in the Enzo convention: a 2D field of total dimensions
(nxt, nyt) is stored flat with index = i + j*nxt, i.e. x (index i) is the
fastest-varying axis.  Fields are Density, TotalEnergy (specific *total* energy
= internal + 0.5 v^2) and Velocity1/2/3; ghosts (NGHOST=3) are zero-gradient
(outflow) on all four sides.

The driver can lay out a Sod shock tube varying along EITHER x or y (uniform in
the transverse direction), so a test can evolve both and assert the resulting
1D profiles are identical -- the rotational-symmetry certification that the y
sweep reproduces the x sweep's physics.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from .. import bridge

NGHOST = 3


def _idx(i, j, nxt):
    return i + j * nxt


def sod_state_2d(axis, nactive, ntrans, gamma, length, disc, left, right):
    """Build the (d, e, u, v, w) flat fields for a Sod tube varying along
    ``axis`` (0 = x, 1 = y), uniform along the transverse axis.

    ``nactive`` active cells along the varying axis, ``ntrans`` along the
    transverse axis.  Returns ``(d, e, u, v, w, dx, nxt, nyt)`` where dx is the
    (uniform, isotropic) cell width and nxt/nyt are the total (active+ghost)
    dimensions.
    """
    dx = length / nactive
    if axis == 0:
        nxa, nya = nactive, ntrans
    else:
        nxa, nya = ntrans, nactive
    nxt = nxa + 2 * NGHOST
    nyt = nya + 2 * NGHOST
    size = nxt * nyt
    d = [0.0] * size
    e = [0.0] * size
    u = [0.0] * size
    v = [0.0] * size
    w = [0.0] * size
    for j in range(nyt):
        for i in range(nxt):
            # active-cell-centred coordinate along the varying axis
            if axis == 0:
                pos = (i - NGHOST + 0.5) * dx
            else:
                pos = (j - NGHOST + 0.5) * dx
            rho, vel, p = left if pos < disc else right
            k = _idx(i, j, nxt)
            d[k] = rho
            # velocity goes into the component for the varying axis
            if axis == 0:
                u[k] = vel
            else:
                v[k] = vel
            e[k] = p / ((gamma - 1.0) * rho) + 0.5 * vel * vel
    return d, e, u, v, w, dx, nxt, nyt


def _fill_ghosts(fields, nxt, nyt):
    """Zero-gradient (outflow) boundaries on all four sides."""
    lo_i, hi_i = NGHOST, nxt - 1 - NGHOST
    lo_j, hi_j = NGHOST, nyt - 1 - NGHOST
    for a in fields:
        # left/right (x) ghosts
        for j in range(nyt):
            for g in range(NGHOST):
                a[_idx(g, j, nxt)] = a[_idx(lo_i, j, nxt)]
                a[_idx(nxt - 1 - g, j, nxt)] = a[_idx(hi_i, j, nxt)]
        # bottom/top (y) ghosts
        for i in range(nxt):
            for g in range(NGHOST):
                a[_idx(i, g, nxt)] = a[_idx(i, lo_j, nxt)]
                a[_idx(i, nyt - 1 - g, nxt)] = a[_idx(i, hi_j, nxt)]


def _max_speed(d, e, u, v, nxt, nyt, gamma):
    cmax = 0.0
    for j in range(NGHOST, nyt - NGHOST):
        for i in range(NGHOST, nxt - NGHOST):
            k = _idx(i, j, nxt)
            vel2 = u[k] * u[k] + v[k] * v[k]
            eint = e[k] - 0.5 * vel2
            p = (gamma - 1.0) * d[k] * eint
            cs = math.sqrt(gamma * max(p, 1e-20) / d[k])
            cmax = max(cmax, math.sqrt(vel2) + cs)
    return cmax


def run(axis=0, t_final=0.2, nactive=128, ntrans=4, cfl=0.4, gamma=1.4,
        length=1.0, disc=0.5, left=(1.0, 0.0, 1.0), right=(0.125, 0.0, 0.1),
        max_steps=100000):
    """Evolve a Sod tube (varying along ``axis``) to ``t_final`` with the full
    2D PPM step.  Returns ``(x, rho, vel, p, nsteps)`` -- a 1D profile sampled
    along the varying axis at the transverse mid-row."""
    d, e, u, v, w, dx, nxt, nyt = sod_state_2d(
        axis, nactive, ntrans, gamma, length, disc, left, right)
    t = 0.0
    n = 0
    dims = [nxt, nyt]
    while t < t_final and n < max_steps:
        _fill_ghosts([d, e, u, v, w], nxt, nyt)
        dt = cfl * dx / _max_speed(d, e, u, v, nxt, nyt, gamma)
        if t + dt > t_final:
            dt = t_final - t
        d, e, u, v, w = bridge.ppm_hydro_step(2, dims, d, e, u, v, w,
                                              dx, dt, gamma)
        t += dt
        n += 1

    # Extract a 1D profile along the varying axis at the transverse midpoint.
    if axis == 0:
        nxa = nactive
        nya = ntrans
        jmid = NGHOST + nya // 2
        x = [(a + 0.5) * dx for a in range(nxa)]
        rho, vel, p = [], [], []
        for a in range(nxa):
            k = _idx(NGHOST + a, jmid, nxt)
            vel2 = u[k] * u[k] + v[k] * v[k]
            eint = e[k] - 0.5 * vel2
            rho.append(d[k]); vel.append(u[k])
            p.append((gamma - 1.0) * d[k] * eint)
    else:
        nya = nactive
        nxa = ntrans
        imid = NGHOST + nxa // 2
        x = [(a + 0.5) * dx for a in range(nya)]
        rho, vel, p = [], [], []
        for a in range(nya):
            k = _idx(imid, NGHOST + a, nxt)
            vel2 = u[k] * u[k] + v[k] * v[k]
            eint = e[k] - 0.5 * vel2
            rho.append(d[k]); vel.append(v[k])
            p.append((gamma - 1.0) * d[k] * eint)
    return x, rho, vel, p, n
