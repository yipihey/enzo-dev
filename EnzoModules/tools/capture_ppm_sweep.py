#!/usr/bin/env python3
"""Capture + validate golden fixtures for the composite PPM 1D sweep.

Two things:

  1. Single-step regression fixtures: take a small shock-tube slice, run one
     ``ppm_sweep_1d``, and store (input slice, params, output slice).  The
     replay test reproduces the output bit-for-bit.

  2. Physics validation (hard gate): evolve a Sod tube to t=0.2 and check the
     L1 density error against the exact Riemann solution.

Usage:
    python3 tools/capture_ppm_sweep.py [--out DIR] [--check-only]
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from enzomodules import bridge, hydro                       # noqa: E402
from enzomodules.fixtures import Fixture, save_fixture      # noqa: E402
from enzomodules.examples import ppm_sod, riemann           # noqa: E402

DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "fixtures", "Hydro", "ppm_sweep_1d")


def single_step_fixture(name, nx, discontinuity, left, right, gamma=1.4,
                        cfl=0.4):
    """Build a deterministic one-sweep fixture from a shock-tube slice."""
    g = ppm_sod.sod_grid(nx=nx, gamma=gamma, discontinuity=discontinuity,
                         left=left, right=right)
    # Boundary fill + pressure, then a deterministic dt (CFL on the IC).
    ppm_sod._fill_ghosts(g)
    p_in = g.pressure()
    dt = cfl * g.dx / ppm_sod.max_wave_speed(g, p_in)

    d0, e0, u0 = list(g.d), list(g.e), list(g.u)
    v0, w0 = list(g.v), list(g.w)
    (d1, e1, u1, v1, w1), (df, ef, uf) = bridge.ppm_sweep_1d(
        d0, e0, u0, v0, w0, p_in,
        g.i1, g.i2, g.dx, dt, gamma, want_fluxes=True)

    fx = Fixture()
    fx.scalars.update(dict(
        kernel="ppm_sweep_1d", name=name, idim=g.idim, i1=g.i1, i2=g.i2,
        nghost=ppm_sod.NGHOST, dx=g.dx, dt=dt, gamma=gamma))
    fx.arrays.update(dict(
        dslice_in=d0, eslice_in=e0, uslice_in=u0, vslice_in=v0, wslice_in=w0,
        pslice_in=p_in,
        dslice_out=d1, eslice_out=e1, uslice_out=u1, vslice_out=v1, wslice_out=w1,
        df=df, ef=ef, uf=uf))
    return fx


SINGLE_STEP_CASES = [
    # name,        nx, disc, left,             right
    ("sod_step",   32, 0.5, (1.0, 0.0, 1.0),  (0.125, 0.0, 0.1)),
    ("strong_step", 32, 0.5, (1.0, 0.0, 10.0), (1.0, 0.0, 0.1)),
]


def validate_sod(nx=200, t_final=0.2, gamma=1.4):
    """Evolve Sod and return the L1 density error vs the exact solution."""
    g, nsteps = ppm_sod.run(t_final=t_final, nx=nx, gamma=gamma)
    x = g.x_centers()
    rho = g.active(g.d)
    rho_ex, _, _ = riemann.sample((1.0, 0.0, 1.0), (0.125, 0.0, 0.1),
                                  gamma, 0.5, x, t_final)
    l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / len(rho)
    return l1, nsteps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()

    bridge.check_precision()
    os.makedirs(args.out, exist_ok=True)

    for name, nx, disc, left, right in SINGLE_STEP_CASES:
        fx = single_step_fixture(name, nx, disc, left, right)
        dmax = max(fx["dslice_out"])
        print(f"single-step {name:12} idim={fx['idim']:3d} dt={fx['dt']:.4e} "
              f"max(rho_out)={dmax:.4f}")
        if not args.check_only:
            save_fixture(os.path.join(args.out, name + ".fixture"), fx)

    l1, nsteps = validate_sod()
    print(f"\nSod evolution to t=0.2: {nsteps} steps, "
          f"L1 density error vs exact = {l1:.4e}")
    # PPM on 200 cells resolves Sod well; this bound catches gross regressions.
    if l1 > 1.5e-2:
        print(f"FAIL: L1 density error {l1:.4e} exceeds 1.5e-2")
        return 1
    print("PPM Sod validation passed.")
    if not args.check_only:
        print(f"Wrote {len(SINGLE_STEP_CASES)} fixtures to {os.path.relpath(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
