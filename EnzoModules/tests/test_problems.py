"""Tests for problem setup / initial-condition generation across Enzo problem
types, via Enzo's InitializeNew (wrapped in enzomodules.problems.Problem).

Requires the grid solver library (deps/build_grid.sh) and the in-repo run/
parameter files.  Skips cleanly otherwise.
"""
import os

import pytest

from enzomodules import problems

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RUN = os.path.join(REPO, "run")


def _param(rel):
    return os.path.join(RUN, rel)


pytestmark = pytest.mark.skipif(
    not problems.available(),
    reason="grid solver library not built; run EnzoModules/deps/build_grid.sh")


# A spread of problem types exercised through the single InitializeNew dispatch.
CASES = [
    ("Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo", 1, 1),
    ("Hydro/Hydro-2D/SedovBlast/SedovBlast.enzo", 7, 2),
    ("Hydro/Hydro-2D/Implosion/Implosion.enzo", 6, 2),
    ("Hydro/Hydro-2D/KelvinHelmholtz/KelvinHelmholtz.enzo", 8, 2),
    ("Hydro/Hydro-2D/NohProblem2D/NohProblem2D.enzo", 9, 2),
]


@pytest.mark.parametrize("rel,expected_type,expected_rank", CASES,
                         ids=[c[0].split("/")[-1] for c in CASES])
def test_problem_initializes(rel, expected_type, expected_rank):
    path = _param(rel)
    if not os.path.exists(path):
        pytest.skip(f"parameter file missing: {rel}")
    with problems.Problem(path) as p:
        assert p.problem_type == expected_type
        assert p.num_grids >= 1
        g = p.grid(0)
        assert g.rank == expected_rank
        assert g.size > 0
        # Every hydro IC carries Density + TotalEnergy + a velocity.
        assert "Density" in g.field_names
        assert "TotalEnergy" in g.field_names
        rho = g.field("Density")
        assert len(rho) == g.size
        assert all(v > 0 for v in rho)        # physical density everywhere


def test_shocktube_initial_state():
    """The Toro-1 (Sod) shock tube IC must have the correct left/right states."""
    path = _param("Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo")
    if not os.path.exists(path):
        pytest.skip("Toro-1 parameter file missing")
    with problems.Problem(path) as p:
        g = p.grid(0)
        nghost = 3
        n = g.size
        rho = g.field("Density")
        active = rho[nghost:n - nghost]
        nx = len(active)
        left = active[nx // 4]
        right = active[3 * nx // 4]
        assert abs(left - 1.0) < 1e-6
        assert abs(right - 0.125) < 1e-6
        # discontinuity: exactly two density levels
        assert abs(max(active) - 1.0) < 1e-6
        assert abs(min(active) - 0.125) < 1e-6


def test_field_access_by_index_and_name():
    path = _param("Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo")
    if not os.path.exists(path):
        pytest.skip("Toro-1 parameter file missing")
    with problems.Problem(path) as p:
        g = p.grid(0)
        by_name = g.field("Density")
        by_index = g.field(g.field_types.index(0))   # Density == FieldType 0
        assert by_name == by_index
