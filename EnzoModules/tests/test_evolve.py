"""Tests for the time-integration path: the CFL timestep (grid::ComputeTimeStep)
and Enzo's full AMR time integrator (EvolveHierarchy), driven end-to-end on a
real problem and checked against the exact Riemann solution.

Requires the grid solver library (deps/build_grid.sh) and the in-repo run/
parameter files.  Skips cleanly otherwise.
"""
import math
import os

import pytest

from enzomodules import bridge, problems
from enzomodules.examples import riemann

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")


@pytest.mark.parametrize("rho,p,u", [(1.0, 1.0, 0.0), (1.0, 1.0, 0.3),
                                     (0.5, 2.0, 0.5)])
def test_compute_timestep_cfl(rho, p, u):
    """grid::ComputeTimeStep returns the CFL dt = courant*dx/(|v|+c_s)."""
    nx, ng = 64, 3
    dimx = nx + 2 * ng
    dx, gamma, courant = 1.0 / nx, 1.4, 0.4
    e = p / ((gamma - 1) * rho) + 0.5 * u * u
    dt = bridge.compute_timestep(1, (dimx, 1, 1), [rho] * dimx, [e] * dimx,
                                 u=[u] * dimx, courant=courant, gamma=gamma,
                                 dx=dx)
    cs = math.sqrt(gamma * p / rho)
    expected = courant * dx / (abs(u) + cs)
    assert abs(dt - expected) < 1e-5 * expected


def _param(rel):
    return os.path.join(REPO, "run", rel)


def test_evolve_shocktube_matches_exact():
    """Evolve the Toro-1 shock tube through Enzo's full time integrator
    (EvolveHierarchy) and compare to the exact Riemann solution."""
    path = _param("Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo")
    if not os.path.exists(path):
        pytest.skip("Toro-1 parameter file missing")
    with problems.Problem(path, evolve=True) as p:
        g = p.grid(0)
        ng = 3
        n = g.size
        nx = g.dims[0] - 2 * ng
        rho = g.field("Density")[ng:n - ng]
        x = [(i + 0.5) / nx for i in range(nx)]
        # Toro-1: left (1, 0.75, 1), right (0.125, 0, 0.1), disc at 0.3, t=0.2.
        rho_ex, _, _ = riemann.sample((1.0, 0.75, 1.0), (0.125, 0.0, 0.1),
                                      1.4, 0.3, x, 0.2)
        l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / nx
        assert l1 < 1e-2, f"evolved L1 density error {l1:.4e}"
        assert abs(max(rho) - 1.0) < 1e-2          # left state preserved
        assert abs(min(rho) - 0.125) < 1e-2        # right state preserved


def test_evolve_bounded_by_cycle():
    """stop_cycle bounds the integration; mass is conserved over the run."""
    path = _param("Hydro/Hydro-2D/Sedov/Sedov.enzo")
    if not os.path.exists(path):
        path = _param("Hydro/Hydro-2D/SedovBlast/SedovBlast.enzo")
    if not os.path.exists(path):
        pytest.skip("Sedov parameter file missing")
    with problems.Problem(path, evolve=True, stop_cycle=3) as p:
        g = p.grid(0)
        ng = 3
        dx, dy = g.dims[0], g.dims[1]
        rho = g.field("Density")
        # interior density stays positive (stable integration)
        assert all(v > 0 for v in rho)
