#!/usr/bin/env python3
"""Capture + validate golden fixtures for the ``twoshock`` kernel.

Runs the legacy Fortran two-shock Riemann solver through the EnzoModules
bridge on a set of canonical Riemann problems, validates the resolved
star-state against an exact solver, and writes golden ``.fixture`` files
consumed by the replay tests.

Usage:
    python3 tools/capture_twoshock.py [--out DIR] [--check-only]
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from enzomodules import bridge                     # noqa: E402
from enzomodules.fixtures import Fixture, save_fixture  # noqa: E402

DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "fixtures", "Hydro", "twoshock")


def run_twoshock(dl, ul, pl, dr, ur, pr, gamma, pmin=1e-20):
    """Single-interface solve; returns (pbar, ubar)."""
    pbar, ubar = bridge.twoshock_raw(
        [dl], [dr], [pl], [pr], [ul], [ur],
        1, 1, 1, 1, 1, 1,
        0.0, gamma, pmin, 0,
        0, [0.0], 0, 0.0,
    )
    return pbar[0], ubar[0]


# ---- exact Riemann solver (Toro Ch.4) for cross-validation -------------

def exact_star(dl, ul, pl, dr, ur, pr, gamma):
    cl = math.sqrt(gamma * pl / dl)
    cr = math.sqrt(gamma * pr / dr)
    g1 = (gamma - 1.0) / (2.0 * gamma)

    def f(p, dk, pk, ck):
        if p > pk:  # shock
            A = 2.0 / ((gamma + 1.0) * dk)
            B = (gamma - 1.0) / (gamma + 1.0) * pk
            return (p - pk) * math.sqrt(A / (p + B))
        return (2.0 * ck) / (gamma - 1.0) * ((p / pk) ** g1 - 1.0)  # rarefaction

    def fp(p):
        return f(p, dl, pl, cl) + f(p, dr, pr, cr) + (ur - ul)

    p = 0.5 * (pl + pr)
    for _ in range(100):
        h = 1e-8 * max(p, 1e-12)
        d = (fp(p + h) - fp(p - h)) / (2 * h)
        pnew = p - fp(p) / d
        if pnew <= 0:
            pnew = 0.5 * p
        if abs(pnew - p) < 1e-12 * pnew:
            p = pnew
            break
        p = pnew
    u = 0.5 * (ul + ur) + 0.5 * (f(p, dr, pr, cr) - f(p, dl, pl, cl))
    return p, u


CASES = [
    # name,            dl,   ul,    pl,  dr,    ur,  pr,   gamma
    ("sod",            1.0,  0.0,   1.0, 0.125, 0.0, 0.1,  1.4),
    ("toro2_123",      1.0, -2.0,   0.4, 1.0,   2.0, 0.4,  1.4),
    ("toro4_collide",  5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950, 1.4),
    ("strong_shock",   1.0,  0.0, 1000.0, 1.0,   0.0, 0.01, 1.4),
    ("symmetric",      2.0,  0.0,   3.0, 2.0,   0.0, 3.0,  1.4),
]


def make_fixture(name, dl, ul, pl, dr, ur, pr, gamma, pbar, ubar):
    fx = Fixture()
    fx.scalars.update(dict(
        kernel="twoshock", name=name, idim=1, jdim=1, i1=1, i2=1, j1=1, j2=1,
        ipresfree=0, gravity=0, idual=0, dt=0.0, gamma=gamma, pmin=1e-20,
        eta1=0.0))
    fx.arrays.update(dict(
        dls=[dl], drs=[dr], pls=[pl], prs=[pr], uls=[ul], urs=[ur],
        grslice=[0.0], pbar=[pbar], ubar=[ubar]))
    return fx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()

    bridge.check_precision()
    os.makedirs(args.out, exist_ok=True)

    print(f"{'case':16} {'pbar(2-shock)':>14} {'pbar(exact)':>12} "
          f"{'ubar(2-shock)':>14} {'ubar(exact)':>12}  note")
    for name, dl, ul, pl, dr, ur, pr, gamma in CASES:
        pbar, ubar = run_twoshock(dl, ul, pl, dr, ur, pr, gamma)
        try:
            pex, uex = exact_star(dl, ul, pl, dr, ur, pr, gamma)
        except Exception:
            pex, uex = float("nan"), float("nan")
        note = ""
        if math.isfinite(pex) and pex > 1e-6:
            relp = abs(pbar - pex) / pex
            note = f"approx (relp={relp:.2f})" if relp > 0.15 else f"ok (relp={relp:.1e})"
        print(f"{name:16} {pbar:14.8g} {pex:12.6g} {ubar:14.8g} {uex:12.6g}  {note}")
        if not args.check_only:
            fx = make_fixture(name, dl, ul, pl, dr, ur, pr, gamma, pbar, ubar)
            save_fixture(os.path.join(args.out, name + ".fixture"), fx)

    # Hard gate: the classic Sod problem (two-shock is accurate here).
    pbar, ubar = run_twoshock(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4)
    pex, uex = exact_star(1.0, 0.0, 1.0, 0.125, 0.0, 0.1, 1.4)
    if abs(pbar - pex) / pex > 0.02 or abs(ubar - uex) > 0.02:
        print(f"\nFAIL: Sod star-state off (pbar={pbar} vs {pex}, ubar={ubar} vs {uex})")
        return 1
    print(f"\nSod star-state validated: pbar={pbar:.6f} (exact {pex:.6f}), "
          f"ubar={ubar:.6f} (exact {uex:.6f})")
    if not args.check_only:
        print(f"Wrote {len(CASES)} fixtures to {os.path.relpath(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
