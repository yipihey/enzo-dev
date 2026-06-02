"""Tests for the wrapped ZEUS solver (grid::ZeusSolver), via the generic
grid-fixture infrastructure.

Requires the grid solver library, which links against the full Enzo shared
library (build with EnzoModules/deps/build_grid.sh).  Skips cleanly otherwise.
"""
import pytest

from enzomodules import bridge
from enzomodules.examples import riemann, zeus_sod

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)


def test_zeus_sod_matches_exact():
    """ZEUS evolves the Sod tube to match the exact Riemann solution.

    ZEUS is more diffusive than the Godunov schemes (artificial viscosity),
    so the tolerance is looser, but density features must still land.
    """
    x, rho, u, p, nsteps = zeus_sod.run(t_final=0.2, nx=200)
    rho_ex, u_ex, p_ex = riemann.sample((1.0, 0.0, 1.0), (0.125, 0.0, 0.1),
                                        1.4, 0.5, x, 0.2)
    n = len(rho)
    l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / n
    assert nsteps > 50
    assert l1 < 1.5e-2, f"ZEUS L1 density error {l1:.3e}"
    assert abs(max(rho) - 1.0) < 1e-2
    assert abs(min(rho) - 0.125) < 1e-2


def test_zeus_single_step_runs():
    """A single ZEUS grid-fixture sweep round-trips fields without error."""
    nghost = 3
    nx = 16
    idim = nx + 2 * nghost
    gamma = 1.4
    dx = 1.0 / nx
    d = [1.0 if (c - nghost) < nx // 2 else 0.125 for c in range(idim)]
    e = [(1.0 if (c - nghost) < nx // 2 else 0.1) / ((gamma - 1.0) * d[c])
         for c in range(idim)]
    u = [0.0] * idim
    d2, e2, u2 = bridge.zeus_sweep_1d(d, e, u, nghost, dx, 1e-3, gamma)
    assert len(d2) == idim
    assert all(v > 0 for v in d2)              # density stays positive
    assert abs(d2[nghost] - 1.0) < 1e-6        # quiescent left edge unchanged
