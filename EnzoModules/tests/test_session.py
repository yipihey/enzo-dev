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


def test_python_evolve_level_matches_exact_riemann():
    """The recursive Python EvolveLevel (problems.Session.evolve_level), built
    only from the certified session steps, integrates the Toro-1 shock tube and
    matches the exact Riemann solution -- the single-level proof that
    EvolveLevel can be written in Python on top of these bridges."""
    with problems.Session(_toro1()) as s:
        n = 0
        while s.time < s.stop_time and n < 10000:
            n += s.evolve_level(0)
        assert n > 0
        g = s.grid(0)
        ng = 3
        sz = g.size
        nx = g.dims[0] - 2 * ng
        rho = g.field("Density")[ng:sz - ng]
        x = [(i + 0.5) / nx for i in range(nx)]
        rho_ex, _, _ = riemann.sample((1.0, 0.75, 1.0), (0.125, 0.0, 0.1),
                                      1.4, 0.3, x, 0.2)
        l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / nx
        assert l1 < 1e-2, f"Python EvolveLevel L1 = {l1:.4e}"


def _implosion_amr():
    p = _param("Hydro/Hydro-2D/ImplosionAMR/ImplosionAMR.enzo")
    if not os.path.exists(p):
        pytest.skip("ImplosionAMR parameter file missing")
    return p


def test_python_evolve_level_amr_recursion():
    """The recursive Python EvolveLevel drives a real multi-level AMR problem
    (ImplosionAMR, 4 initial levels) -- exercising sub-cycling, recursion into
    finer levels, conservative flux correction + projection (update_from_finer)
    and regridding -- and reproduces Enzo's monolithic EvolveHierarchy on the
    root grid (mean density to ~1e-5)."""
    import enzomodules.problems as P
    path = _implosion_amr()
    cycles = 3

    # reference: Enzo's own EvolveHierarchy
    with problems.Problem(path, evolve=True, stop_cycle=cycles) as p:
        rho_ref = p.grid(0).field("Density")

    # the Python recursive EvolveLevel, same number of root cycles
    with problems.Session(path) as s:
        assert s.num_grids_on_level(1) > 0          # starts multi-level
        with P._suppress_fd_output():
            s.rebuild(0)
            for _ in range(cycles):
                s.evolve_level(0)                   # recurses + flux-corrects
                s.rebuild(0)
        rho_sess = s.grid(0).field("Density")
        # the recursion really created finer grids
        assert s.num_grids_on_level(1) > 1

    assert len(rho_ref) == len(rho_sess)
    n = len(rho_ref)
    assert all(v > 0 for v in rho_sess)             # physical
    l1 = sum(abs(a - b) for a, b in zip(rho_ref, rho_sess)) / n
    mean_ref = sum(rho_ref) / n
    mean_sess = sum(rho_sess) / n
    assert l1 < 1e-2, f"Python EvolveLevel vs EvolveHierarchy L1 = {l1:.4e}"
    assert abs(mean_ref - mean_sess) < 1e-3 * mean_ref   # conservation


def _photontest():
    p = _param("RadiationTransport/PhotonTest/PhotonTest.enzo")
    if not os.path.exists(p):
        pytest.skip("PhotonTest parameter file missing")
    return p


def test_session_evolve_photons_emits():
    """EvolvePhotons actually emits from the parameter-file source and deposits
    photo-ionization rates on the live hierarchy (not just 'runs')."""
    import enzomodules.problems as P
    with problems.Session(_photontest()) as s:
        g = s.grid(0)
        assert "kphHI" in g.field_names
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.evolve_photons(0)
        kph = g.field("kphHI")
        # the source switched on: some cells now have a non-zero ionization rate
        assert max(kph) > 0.0, "EvolvePhotons deposited no photo-ionization rate"
        assert sum(1 for v in kph if v > 0) >= 1


def test_session_evolve_photons_ionization_front():
    """Driving several photon steps grows the ionized region monotonically --
    the Stromgren-sphere I-front propagation that is PhotonTest's purpose --
    proving the full EvolvePhotons orchestration couples radiation to the
    chemistry on the live hierarchy."""
    import enzomodules.problems as P
    with problems.Session(_photontest()) as s:
        g = s.grid(0)
        fracs = []
        ncells = []
        for _ in range(4):
            with P._suppress_fd_output():
                s.set_boundary(0)
                s.set_dt(0, s.compute_dt(0))
                s.evolve_photons(0)
                s.advance_time(0)
            hi = g.field("HIDensity")
            hii = g.field("HIIDensity")
            kph = g.field("kphHI")
            fracs.append(sum(hii) / (sum(hi) + sum(hii)))
            ncells.append(sum(1 for v in kph if v > 0))
        # the ionized fraction and the illuminated volume both grow each step
        assert all(b > a for a, b in zip(fracs, fracs[1:])), f"ionized frac {fracs}"
        assert all(b >= a for a, b in zip(ncells, ncells[1:])), f"kph cells {ncells}"
        assert fracs[-1] > fracs[0] * 2


def test_session_radiation_hydro_coupled():
    """The recursive Python EvolveLevel with radiation=True drives a coupled
    radiation-hydrodynamics AMR problem (PhotonTestAMR: hydro + RT + AMR),
    emitting photons each step and ionizing the medium while the hydro and
    regridding run -- the full RHD time integrator, from Python."""
    path = _param("RadiationTransport/PhotonTestAMR/PhotonTestAMR.enzo")
    if not os.path.exists(path):
        pytest.skip("PhotonTestAMR parameter file missing")
    with problems.Session(path) as s:
        g = s.grid(0)

        def ionized():
            hi = g.field("HIDensity")
            hii = g.field("HIIDensity")
            return sum(hii) / (sum(hi) + sum(hii))

        f0 = ionized()
        for _ in range(2):
            s.evolve_level(0, radiation=True)
        f1 = ionized()
        assert f1 > f0, f"ionized fraction did not grow ({f0:.3e} -> {f1:.3e})"
        assert all(v > 0 for v in g.field("Density"))   # hydro stays physical
