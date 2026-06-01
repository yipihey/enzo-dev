#!/usr/bin/env python3
"""Run the Sod shock tube through the wrapped legacy PPM solver and compare
to the exact Riemann solution -- a self-contained demonstration that the
extracted kernel works end-to-end.

    python3 tools/demo_sod.py [--nx 200] [--t 0.2]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from enzomodules.examples import ppm_sod, riemann          # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=200)
    ap.add_argument("--t", type=float, default=0.2)
    args = ap.parse_args()

    gamma = 1.4
    g, nsteps = ppm_sod.run(t_final=args.t, nx=args.nx, gamma=gamma)
    x = g.x_centers()
    rho, u, p = g.active(g.d), g.active(g.u), g.active(g.pressure())
    rho_e, u_e, p_e = riemann.sample((1.0, 0.0, 1.0), (0.125, 0.0, 0.1),
                                     gamma, 0.5, x, args.t)

    n = len(rho)
    l1 = sum(abs(a - b) for a, b in zip(rho, rho_e)) / n
    print(f"Sod shock tube: nx={args.nx}, t={args.t}, steps={nsteps}, "
          f"L1(rho)={l1:.3e}\n")
    print(f"{'x':>6} {'rho':>8} {'rho_exact':>10} {'u':>8} {'p':>8}   density (PPM=#, exact=.)")
    rmin, rmax = 0.0, 1.05
    for k in range(0, n, max(1, n // 32)):
        col = int((rho[k] - rmin) / (rmax - rmin) * 40)
        cole = int((rho_e[k] - rmin) / (rmax - rmin) * 40)
        bar = ["."] * 42
        bar[min(41, max(0, cole))] = "."
        bar[min(41, max(0, col))] = "#"
        print(f"{x[k]:6.3f} {rho[k]:8.4f} {rho_e[k]:10.4f} {u[k]:8.4f} "
              f"{p[k]:8.4f}   {''.join(bar)}")


if __name__ == "__main__":
    main()
