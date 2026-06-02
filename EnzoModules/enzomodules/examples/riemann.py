"""Exact Riemann solver for the 1D Euler equations (Toro, Ch. 4).

Independent of Enzo -- this is the analytic *truth* the wrapped PPM kernels
are validated against.  Handles shocks and rarefactions on both sides and can
sample the full self-similar solution W(x, t), so tests can compute an L1
error over the whole domain rather than eyeballing a few points.
"""
from __future__ import annotations

import math
from typing import Callable, Tuple

State = Tuple[float, float, float]  # (density, velocity, pressure)


def _f_and_df(p: float, dk: float, pk: float, ck: float, gamma: float):
    """Pressure function for one side and its derivative (Toro eq. 4.6/4.7)."""
    if p > pk:  # shock
        A = 2.0 / ((gamma + 1.0) * dk)
        B = (gamma - 1.0) / (gamma + 1.0) * pk
        q = math.sqrt(A / (p + B))
        f = (p - pk) * q
        df = q * (1.0 - 0.5 * (p - pk) / (p + B))
    else:       # rarefaction
        g1 = (gamma - 1.0) / (2.0 * gamma)
        f = (2.0 * ck) / (gamma - 1.0) * ((p / pk) ** g1 - 1.0)
        df = (1.0 / (dk * ck)) * (p / pk) ** (-(gamma + 1.0) / (2.0 * gamma))
    return f, df


def star_state(WL: State, WR: State, gamma: float) -> Tuple[float, float]:
    """Solve for the star-region pressure and velocity (p*, u*)."""
    dL, uL, pL = WL
    dR, uR, pR = WR
    cL = math.sqrt(gamma * pL / dL)
    cR = math.sqrt(gamma * pR / dR)

    p = 0.5 * (pL + pR)
    for _ in range(100):
        fL, dfL = _f_and_df(p, dL, pL, cL, gamma)
        fR, dfR = _f_and_df(p, dR, pR, cR, gamma)
        f = fL + fR + (uR - uL)
        p_new = p - f / (dfL + dfR)
        if p_new <= 0:
            p_new = 0.5 * p
        if abs(p_new - p) < 1e-12 * (p_new + p):
            p = p_new
            break
        p = p_new
    fL, _ = _f_and_df(p, dL, pL, cL, gamma)
    fR, _ = _f_and_df(p, dR, pR, cR, gamma)
    u = 0.5 * (uL + uR) + 0.5 * (fR - fL)
    return p, u


def solver(WL: State, WR: State, gamma: float) -> Callable[[float], State]:
    """Return ``W(s)`` giving the exact state at self-similar coordinate
    ``s = x / t`` (with the discontinuity at ``x = 0``)."""
    dL, uL, pL = WL
    dR, uR, pR = WR
    cL = math.sqrt(gamma * pL / dL)
    cR = math.sqrt(gamma * pR / dR)
    pstar, ustar = star_state(WL, WR, gamma)
    g = gamma

    def W(s: float) -> State:
        if s <= ustar:  # left of contact
            if pstar > pL:  # left shock
                qL = math.sqrt((pstar / pL + (g - 1) / (g + 1)) / (2 * g / (g + 1)))
                SL = uL - cL * qL
                if s <= SL:
                    return WL
                dstar = dL * ((pstar / pL + (g - 1) / (g + 1)) /
                              ((g - 1) / (g + 1) * pstar / pL + 1))
                return (dstar, ustar, pstar)
            else:           # left rarefaction
                cstarL = cL * (pstar / pL) ** ((g - 1) / (2 * g))
                SHL = uL - cL
                STL = ustar - cstarL
                if s <= SHL:
                    return WL
                if s >= STL:
                    dstarL = dL * (pstar / pL) ** (1 / g)
                    return (dstarL, ustar, pstar)
                # inside fan
                u = 2 / (g + 1) * (cL + (g - 1) / 2 * uL + s)
                c = 2 / (g + 1) * (cL + (g - 1) / 2 * (uL - s))
                d = dL * (c / cL) ** (2 / (g - 1))
                p = pL * (c / cL) ** (2 * g / (g - 1))
                return (d, u, p)
        else:           # right of contact
            if pstar > pR:  # right shock
                qR = math.sqrt((pstar / pR + (g - 1) / (g + 1)) / (2 * g / (g + 1)))
                SR = uR + cR * qR
                if s >= SR:
                    return WR
                dstar = dR * ((pstar / pR + (g - 1) / (g + 1)) /
                              ((g - 1) / (g + 1) * pstar / pR + 1))
                return (dstar, ustar, pstar)
            else:           # right rarefaction
                cstarR = cR * (pstar / pR) ** ((g - 1) / (2 * g))
                SHR = uR + cR
                STR = ustar + cstarR
                if s >= SHR:
                    return WR
                if s <= STR:
                    dstarR = dR * (pstar / pR) ** (1 / g)
                    return (dstarR, ustar, pstar)
                u = 2 / (g + 1) * (-cR + (g - 1) / 2 * uR + s)
                c = 2 / (g + 1) * (cR - (g - 1) / 2 * (uR - s))
                d = dR * (c / cR) ** (2 / (g - 1))
                p = pR * (c / cR) ** (2 * g / (g - 1))
                return (d, u, p)

    return W


def sample(WL: State, WR: State, gamma: float, x0: float, x, t: float):
    """Exact ``(rho, u, p)`` arrays at positions ``x`` (iterable) and time
    ``t`` for an initial discontinuity at ``x0``."""
    W = solver(WL, WR, gamma)
    rho, u, p = [], [], []
    for xi in x:
        d, v, pr = W((xi - x0) / t) if t > 0 else (WL if xi < x0 else WR)
        rho.append(d); u.append(v); p.append(pr)
    return rho, u, p
