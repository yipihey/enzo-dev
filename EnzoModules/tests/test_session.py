"""Tests for the persistent-hierarchy *session*: a live Enzo hierarchy whose
time loop is driven step-by-step from Python (problems.Session), rather than
handed wholesale to EvolveHierarchy.

The headline guarantee is equivalence: a Python-driven loop calling the
certified legacy steps (SetBoundaryConditions -> ComputeTimeStep ->
SolveHydroEquations -> advance) reproduces Enzo's monolithic EvolveHierarchy
*bit-for-bit* on a shock tube, and matches the exact Riemann solution.  The
gravity / particle / radiation steps are exercised for "runs on the live
hierarchy"; their physics is certified in the dedicated bridge tests.

Requires the grid solver library (deps/build_grid.sh) and the in-repo run/
parameter files.  Skips cleanly otherwise.
"""
import os

import pytest

from enzomodules import bridge, problems
from enzomodules.examples import riemann

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")


def _param(rel):
    return os.path.join(REPO, "run", rel)


def _toro1():
    p = _param("Hydro/Hydro-1D/Toro-1-ShockTube/Toro-1-ShockTube.enzo")
    if not os.path.exists(p):
        pytest.skip("Toro-1 parameter file missing")
    return p


def test_session_matches_exact_riemann():
    """A Python-driven timestep loop on the Toro-1 shock tube reproduces the
    exact Riemann solution (same L1 as EvolveHierarchy)."""
    with problems.Session(_toro1()) as s:
        assert s.num_grids == 1
        assert s.time == 0.0
        cycles = s.run(level=0)
        assert cycles > 0
        assert abs(s.time - s.stop_time) < 1e-6 * s.stop_time

        g = s.grid(0)
        ng = 3
        n = g.size
        nx = g.dims[0] - 2 * ng
        rho = g.field("Density")[ng:n - ng]
        x = [(i + 0.5) / nx for i in range(nx)]
        rho_ex, _, _ = riemann.sample((1.0, 0.75, 1.0), (0.125, 0.0, 0.1),
                                      1.4, 0.3, x, 0.2)
        l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / nx
        assert l1 < 1e-2, f"session L1 density error {l1:.4e}"


def test_session_reproduces_evolvehierarchy_bitwise():
    """The Python-driven loop and Enzo's monolithic EvolveHierarchy produce
    identical state -- the equivalence that lets EvolveLevel be re-implemented
    on top of these certified steps."""
    path = _toro1()
    with problems.Problem(path, evolve=True) as p:
        rho_mono = p.grid(0).field("Density")
    with problems.Session(path) as s:
        s.run(level=0)
        rho_sess = s.grid(0).field("Density")
    assert len(rho_mono) == len(rho_sess)
    linf = max(abs(a - b) for a, b in zip(rho_mono, rho_sess))
    assert linf == 0.0, f"session vs EvolveHierarchy Linf = {linf:.3e}"


def test_session_manual_steps():
    """The individual orchestration steps can be called by hand, and a manual
    loop advances time to StopTime with positive density throughout."""
    with problems.Session(_toro1()) as s:
        n = 0
        while s.time < s.stop_time and n < 10000:
            s.set_boundary(0)
            dt = s.compute_dt(0)
            assert dt > 0
            s.set_dt(0, dt)
            s.solve_hydro(0)
            s.advance_time(0)
            n += 1
        assert n > 0
        rho = s.grid(0).field("Density")
        assert all(v > 0 for v in rho)


def test_session_gravity_chain_runs():
    """The self-gravity chain (deposit + Poisson solve + accelerations) runs on
    the live hierarchy of a self-gravitating problem."""
    path = _param("GravitySolver/GravityTest/GravityTest.enzo")
    if not os.path.exists(path):
        pytest.skip("GravityTest parameter file missing")
    with problems.Session(path) as s:
        s.set_boundary(0)
        s.gravity(0)          # raises on FAIL
        # particle drift is callable on the live hierarchy
        if s.grid(0).num_particles > 0:
            dt = s.compute_dt(0)
            s.set_dt(0, dt)
            s.update_particles(0)


def test_session_evolve_photons_runs():
    """The radiative-transfer step is callable on a live RT hierarchy (the
    transport physics itself is certified in test_radiation.py)."""
    path = _param("RadiationTransport/PhotonTest/PhotonTest.enzo")
    if not os.path.exists(path):
        pytest.skip("PhotonTest parameter file missing")
    with problems.Session(path) as s:
        s.set_boundary(0)
        dt = s.compute_dt(0)
        s.set_dt(0, dt)
        s.evolve_photons(0)   # raises on FAIL
        # RadiativeTransferInitialize allocated the RT fields
        assert "kphHI" in s.grid(0).field_names
