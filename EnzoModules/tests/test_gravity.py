"""Tests for the gravity/Poisson solver (Enzo's MultigridSolver, the engine
behind grid::SolveForPotential).

Validated with a manufactured solution: phi = product of sin(pi x_d), which
vanishes on the (Dirichlet) boundary and has continuous Laplacian
-(rank) pi^2 phi.  Solving with that rhs must recover phi's shape (the solver
carries an internal cell-size normalization, so we compare shape via a fitted
scale).  Requires the grid solver library; skips otherwise.
"""
import math

import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)


def _manufactured(dims):
    """phi_exact and rhs = continuous-Laplacian(phi) for phi = prod sin(pi x_d)."""
    rank = sum(1 for v in dims if v > 1)
    ns = [d if d > 1 else 1 for d in dims]
    nx, ny, nz = (ns + [1, 1])[:3]
    phi = [0.0] * (nx * ny * nz)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = (k * ny + j) * nx + i
                v = math.sin(math.pi * i / (nx - 1))
                if ny > 1:
                    v *= math.sin(math.pi * j / (ny - 1))
                if nz > 1:
                    v *= math.sin(math.pi * k / (nz - 1))
                phi[idx] = v
    rhs = [-rank * math.pi ** 2 * v for v in phi]
    return phi, rhs


def _shape_relerr(sol, phi):
    interior = [i for i, v in enumerate(phi) if abs(v) > 1e-6]
    num = sum(sol[i] * phi[i] for i in interior)
    den = sum(phi[i] ** 2 for i in interior)
    scale = num / den
    peak = max(abs(scale * phi[i]) for i in interior)
    return max(abs(sol[i] - scale * phi[i]) for i in interior) / peak, scale


@pytest.mark.parametrize("dims", [(33, 1, 1), (33, 33, 1), (17, 17, 17)])
def test_manufactured_solution(dims):
    phi, rhs = _manufactured(dims)
    sol, norm, mean = bridge.poisson_solve(rhs, dims)
    # Solver converged.
    assert mean > 0 and norm / mean < 1e-5
    # Recovered the manufactured solution's shape (4th-order accurate).
    relerr, scale = _shape_relerr(sol, phi)
    assert relerr < 5e-3, f"shape relerr {relerr:.2e} for dims {dims}"
    assert abs(scale) > 0


def test_linearity():
    dims = (33, 1, 1)
    _, rhs = _manufactured(dims)
    s1, _, _ = bridge.poisson_solve(rhs, dims)
    s2, _, _ = bridge.poisson_solve([2.0 * r for r in rhs], dims)
    for a, b in zip(s1, s2):
        if abs(a) > 1e-8:
            assert abs(b / a - 2.0) < 1e-4


def test_superposition():
    dims = (33, 1, 1)
    n = dims[0]
    rhs1 = [math.sin(math.pi * i / (n - 1)) for i in range(n)]
    rhs2 = [math.sin(2 * math.pi * i / (n - 1)) for i in range(n)]
    s1, _, _ = bridge.poisson_solve(rhs1, dims)
    s2, _, _ = bridge.poisson_solve(rhs2, dims)
    s12, _, _ = bridge.poisson_solve([a + b for a, b in zip(rhs1, rhs2)], dims)
    for a, b, c in zip(s1, s2, s12):
        assert abs(c - (a + b)) < 1e-4 * (abs(a) + abs(b) + 1e-6)
