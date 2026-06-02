"""Tests for the wrapped full multi-dimensional PPM hydro step
(grid::SolveHydroEquations, which runs the x/y/z directional Euler sweeps),
via the generic grid-fixture infrastructure.

The key certification that the y (and, by symmetry, z) sweeps work is
ROTATIONAL SYMMETRY: a Sod tube varying along x and the SAME tube varying
along y, evolved identically, must produce bitwise-near-identical 1D profiles.
We also compare the evolved density profile to the exact Riemann solution.

Requires the grid solver library, which links against the full Enzo shared
library (build with EnzoModules/deps/build_grid.sh).  Skips cleanly otherwise.
"""
import pytest

from enzomodules import bridge
from enzomodules.examples import ppm_hydro2d, riemann

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)

WL = (1.0, 0.0, 1.0)
WR = (0.125, 0.0, 0.1)
GAMMA = 1.4
T_FINAL = 0.2
NACTIVE = 128


def test_ppm_xy_rotational_symmetry():
    """The x sweep and y sweep must produce identical physics.

    Evolve the same Sod tube along x and along y; the 1D profiles must match
    to ~1e-10 -- this certifies grid::yEulerSweep against grid::xEulerSweep.
    """
    xx, rho_x, u_x, p_x, nx = ppm_hydro2d.run(
        axis=0, t_final=T_FINAL, nactive=NACTIVE, left=WL, right=WR, gamma=GAMMA)
    xy, rho_y, u_y, p_y, ny = ppm_hydro2d.run(
        axis=1, t_final=T_FINAL, nactive=NACTIVE, left=WL, right=WR, gamma=GAMMA)

    assert nx == ny and nx > 20
    drho = max(abs(a - b) for a, b in zip(rho_x, rho_y))
    du = max(abs(a - b) for a, b in zip(u_x, u_y))
    dp = max(abs(a - b) for a, b in zip(p_x, p_y))
    assert drho < 1e-10, f"x/y density profiles differ by {drho:.3e}"
    assert du < 1e-10, f"x/y velocity profiles differ by {du:.3e}"
    assert dp < 1e-10, f"x/y pressure profiles differ by {dp:.3e}"


def test_ppm_sod_matches_exact():
    """The full PPM step evolves the Sod tube to match the exact Riemann
    solution (small L1 density error on a resolved grid)."""
    x, rho, u, p, nsteps = ppm_hydro2d.run(
        axis=0, t_final=T_FINAL, nactive=NACTIVE, left=WL, right=WR, gamma=GAMMA)
    rho_ex, u_ex, p_ex = riemann.sample(WL, WR, GAMMA, 0.5, x, T_FINAL)
    n = len(rho)
    l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / n
    assert nsteps > 20
    assert l1 < 2e-2, f"PPM L1 density error {l1:.3e}"
    assert abs(max(rho) - 1.0) < 2e-2
    assert abs(min(rho) - 0.125) < 2e-2
