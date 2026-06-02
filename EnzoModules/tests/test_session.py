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


def _interior(g, field):
    """Field values with the 3-zone ghost boundary stripped (1D helper)."""
    ng, n = 3, g.dims[0]
    return list(g.field(field))[ng:n - ng]


def test_session_checkpoint_restart_bitwise(tmp_path):
    """write_output -> from_output is a true checkpoint/restart: running N steps,
    checkpointing, reloading and continuing M steps reproduces an uninterrupted
    N+M-step run *bit-for-bit* (the gold-standard restart guarantee)."""
    import enzomodules.problems as P
    # uninterrupted reference: 6 steps
    with problems.Session(_toro1()) as ref:
        with P._suppress_fd_output():
            for _ in range(6):
                ref.step(0)
        rho_ref = _interior(ref.grid(0), "Density")
        te_ref = _interior(ref.grid(0), "TotalEnergy")
        t_ref, cyc_ref = ref.time, ref.cycle
    # 3 steps -> checkpoint (to a persistent dir) -> reload -> 3 more steps
    with problems.Session(_toro1(), workdir=str(tmp_path)) as s:
        with P._suppress_fd_output():
            for _ in range(3):
                s.step(0)
        path = s.write_output(number=0)
    s2 = problems.Session.from_output(path)
    try:
        with P._suppress_fd_output():
            for _ in range(3):
                s2.step(0)
        assert s2.cycle == cyc_ref
        assert abs(s2.time - t_ref) <= 1e-12 * max(1.0, abs(t_ref))
        assert _interior(s2.grid(0), "Density") == rho_ref      # bitwise
        assert _interior(s2.grid(0), "TotalEnergy") == te_ref
    finally:
        s2.close()


def test_session_write_output_multigrid(tmp_path):
    """write_output / from_output preserve a multi-grid AMR hierarchy: the grid
    count and the root-grid interior state round-trip exactly."""
    path = _param("Hydro/Hydro-2D/ImplosionAMR/ImplosionAMR.enzo")
    if not os.path.exists(path):
        pytest.skip("ImplosionAMR parameter file missing")
    import enzomodules.problems as P
    with problems.Session(path, workdir=str(tmp_path)) as s:
        s.run_amr(max_cycles=3)              # run_amr self-suppresses output
        with P._suppress_fd_output():
            s.set_boundary(0)                # consistent ghost zones for compare
        ngrids = s.num_grids
        rho = list(s.grid(0).field("Density"))
        dump = s.write_output(number=0)
    s2 = problems.Session.from_output(dump)
    try:
        assert s2.num_grids == ngrids                # hierarchy preserved
        with P._suppress_fd_output():
            s2.set_boundary(0)
        assert list(s2.grid(0).field("Density")) == rho   # bitwise roundtrip
    finally:
        s2.close()


def test_session_from_output_failure_is_graceful():
    """Hardening: from_output on a nonexistent dump raises instead of aborting."""
    import tempfile
    with pytest.raises(RuntimeError):
        problems.Session.from_output(
            os.path.join(tempfile.gettempdir(), "enzomodules_no_such_dump0000"))


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


def test_session_mhd_rk_runs():
    """solve_hydro dispatches the Runge-Kutta MHD path (HD_RK / MHD_RK) on the
    live hierarchy: the 2nd-order two-step integration (1st step, refresh
    boundaries, 2nd step) with Dedner wave speeds and NColor set up first.  Run
    the Brio-Wu MHD shock tube a few dozen cycles and check the solution stays
    finite with the density (mass) conserved -- the RK MHD solver driven
    step-by-step from Python, not just the kernel line-solvers.  (HD_RK uses the
    same two-step machinery with RungeKutta2 in place of MHDRK2.)"""
    path = _param("MHD/1D/BrioWu-MHD-1D/BrioWu-MHD-1D.enzo")
    if not os.path.exists(path):
        pytest.skip("BrioWu-MHD-1D parameter file missing")
    with problems.Session(path) as s:
        mass0 = sum(s.grid(0).field("Density"))
        n = s.run_amr(max_cycles=40)
        rho = s.grid(0).field("Density")
        assert n > 0
        assert all(v == v and abs(v) < 1e30 for v in rho)   # finite, no NaN/Inf
        assert abs(sum(rho) - mass0) < 1e-6 * mass0          # mass conserved


def test_session_cosmology_expansion():
    """The cosmology accessor (CosmologyComputeExpansionFactor) reports the scale
    factor / redshift on a comoving run, and (1, 0) on a non-comoving one."""
    path = _param("Cosmology/SphericalInfall/SphericalInfall.enzo")
    if not os.path.exists(path):
        pytest.skip("SphericalInfall parameter file missing")
    with problems.Session(path) as s:
        a, z = s.cosmology()
        # Enzo normalizes a = 1 at the initial redshift (here z = 99).
        assert abs(a - 1.0) < 1e-6
        assert abs(z - 99.0) < 1e-3
        assert s.scale_factor == a and s.redshift == z
    # non-cosmological problem -> (1, 0)
    with problems.Session(_toro1()) as s:
        assert s.cosmology() == (1.0, 0.0)


def test_session_physics_steps_noop():
    """Hardening: the optional per-step physics modules (active particles, UV
    background, turbulence forcing, conduction, shock finding) are clean no-ops
    on a plain hydro problem with that physics off, leaving the state intact."""
    import enzomodules.problems as P
    with problems.Session(_toro1()) as s:
        rho0 = list(s.grid(0).field("Density"))
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.active_particles(0)
            s.update_radiation_field(0)
            s.random_forcing(0)
            s.conduct_heat(0)
            s.find_shocks(0)
        assert list(s.grid(0).field("Density")) == rho0


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


def test_session_inline_halo_finder():
    """Enzo's inline FOF halo finder runs on the live particle hierarchy
    (GravityTest, 5000 particles) and returns a consistent catalogue: at least
    one halo, each at least min_size particles, all grouped particles within the
    total, centre of mass inside the box."""
    path = _param("GravitySolver/GravityTest/GravityTest.enzo")
    if not os.path.exists(path):
        pytest.skip("GravityTest parameter file missing")
    import enzomodules.problems as P
    with problems.Session(path) as s:
        with P._suppress_fd_output():
            halos = s.find_halos(linking_length=0.2, min_size=32)
        total = sum(s.grid(i).num_particles for i in range(s.num_grids))
        assert total == 5000
        assert len(halos) >= 1
        assert halos == sorted(halos, key=lambda h: -h["n"])   # largest first
        grouped = sum(h["n"] for h in halos)
        assert 0 < grouped <= total
        for h in halos:
            assert h["n"] >= 32
            assert all(0.0 <= c <= 1.0 for c in h["center"])
        # a larger linking length percolates more particles into halos
        with P._suppress_fd_output():
            halos_b = s.find_halos(linking_length=0.6, min_size=32)
        assert sum(h["n"] for h in halos_b) >= grouped


def test_session_inline_halo_finder_idempotent():
    """The halo finder moves particles off the grids and restores them, so it is
    repeatable: particle count is conserved, the catalogue is stable across
    calls, and the session can still be evolved afterwards."""
    path = _param("GravitySolver/GravityTest/GravityTest.enzo")
    if not os.path.exists(path):
        pytest.skip("GravityTest parameter file missing")
    import enzomodules.problems as P
    with problems.Session(path) as s:
        counts = []
        with P._suppress_fd_output():
            for _ in range(3):
                halos = s.find_halos(linking_length=0.2, min_size=32)
                counts.append((len(halos), sum(s.grid(i).num_particles
                                               for i in range(s.num_grids))))
        assert all(c == counts[0] for c in counts)             # stable
        assert counts[0][1] == 5000                            # conserved
        # the session is still usable after halo finding
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.solve_hydro(0)
            s.advance_time(0)
        assert s.time > 0.0


def test_session_subfind():
    """SUBFIND finds subgroups inside the FOF halos (the big GravityTest clump
    has >= 2*DesLinkNgb particles), reported via subhalo_count."""
    path = _param("GravitySolver/GravityTest/GravityTest.enzo")
    if not os.path.exists(path):
        pytest.skip("GravityTest parameter file missing")
    import enzomodules.problems as P
    with problems.Session(path) as s:
        with P._suppress_fd_output():
            halos = s.find_halos(subfind=True, linking_length=0.2, min_size=32)
        assert len(halos) >= 1
        assert s.subhalo_count() >= 1          # at least one resolved subgroup


def test_session_halo_finder_no_particles():
    """Hardening: the halo finder is a clean no-op (returns []) on a problem
    with no particles, and does not disturb the hydro state."""
    import enzomodules.problems as P
    with problems.Session(_toro1()) as s:
        assert s.grid(0).num_particles == 0
        with P._suppress_fd_output():
            assert s.find_halos() == []
            assert s.find_halos(subfind=True) == []
            s.run(0)                            # hydro unaffected
        assert all(v > 0 for v in s.grid(0).field("Density"))


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


def test_session_evolve_photons_star_sources():
    """The star-particle radiation path (StarParticleInitialize ->
    StarParticleRadTransfer -> EvolvePhotons) runs on a live hierarchy with a
    star particle (ProblemType 252).  This used to segfault (no SubgridMarker,
    unsized dtPhoton); it now completes cleanly.  Whether the star emits depends
    on Enzo's formation lifecycle marking it an active radiation source, so we
    assert the pipeline is robust rather than a specific deposited rate."""
    import enzomodules.problems as P
    path = _param("StarParticle/RadiatingStarParticleSingleTest/"
                  "TestRadiatingStarParticleSingle.enzo")
    if not os.path.exists(path):
        pytest.skip("RadiatingStarParticleSingleTest parameter file missing")
    with problems.Session(path) as s:
        g = s.grid(0)
        assert g.num_particles >= 1                  # the star particle exists
        assert "kphHI" in g.field_names              # RT fields allocated
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.evolve_photons(0, stars=True)          # crash-free (raises on FAIL)
            kph = g.field("kphHI")
        assert all(v >= 0 for v in kph)              # rates stay non-negative


def test_session_star_particle_lifecycle_radiates():
    """The full star-particle lifecycle makes a star form, activate and RADIATE.
    Driving ProblemType 252 (a single Pop III star) with star_particles
    (StarParticleInitialize -> StarParticleHandler -> StarParticleFinalize /
    ActivateNewStar) plus evolve_photons(stars=True), at timesteps finer than the
    star's lifetime, the star activates (type flips from unborn to active) and
    ionizes the surrounding gas: a Stromgren I-front whose illuminated volume
    grows monotonically.  This is the end-to-end proof that the formation ->
    feedback -> activation -> radiation chain works (previously the star never
    activated, so evolve_photons(stars=True) was a no-op)."""
    import enzomodules.problems as P
    path = _param("StarParticle/RadiatingStarParticleSingleTest/"
                  "TestRadiatingStarParticleSingle.enzo")
    if not os.path.exists(path):
        pytest.skip("RadiatingStarParticleSingleTest parameter file missing")
    with problems.Session(path) as s:
        g = s.grid(0)
        ncells = []
        with P._suppress_fd_output():
            for _ in range(8):
                s.set_boundary(0)
                # timestep finer than the star's (short Pop III) lifetime, so we
                # resolve its radiating main-sequence window
                s.set_dt(0, min(s.compute_dt(0), 0.03))
                s.star_particles(0)                  # form / feedback / activate
                s.evolve_photons(0, stars=True)      # radiate from active star
                s.advance_time(0)
                kph = g.field("kphHI")
                ncells.append(sum(1 for v in kph if v > 0))
        assert max(ncells) > 0                       # the star radiated
        assert ncells[-1] > ncells[0]                # I-front propagated outward
        assert ncells == sorted(ncells)              # illuminated volume grew monotonically


def test_session_star_particles_noop():
    """Hardening: star_particles is a clean no-op on a problem with no star
    physics (Toro-1), leaving the state intact."""
    import enzomodules.problems as P
    with problems.Session(_toro1()) as s:
        rho0 = list(s.grid(0).field("Density"))
        with P._suppress_fd_output():
            s.set_dt(0, s.compute_dt(0))
            s.star_particles(0)                      # no-op
        assert list(s.grid(0).field("Density")) == rho0


def test_session_evolve_photons_noop_without_rt():
    """Hardening: evolve_photons (both source modes) is a clean no-op on a
    problem with RadiativeTransfer off, and leaves the hydro path intact."""
    import enzomodules.problems as P
    with problems.Session(_toro1()) as s:
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.evolve_photons(0)                      # no RT -> no-op, no crash
            s.evolve_photons(0, stars=True)          # no-op, no crash
            s.run(0)                                 # hydro unaffected
        assert abs(s.time - s.stop_time) < 1e-6 * s.stop_time
        assert all(v > 0 for v in s.grid(0).field("Density"))


def test_session_solve_cooling():
    """The radiative-cooling + chemistry sub-step (grid::MultiSpeciesHandler)
    runs on the live hierarchy: after heating/ionizing the gas with a few RT
    steps, repeated solve_cooling (no further heating) monotonically lowers the
    total energy and evolves the species -- real, certified cooling/chemistry,
    not just a callable no-op."""
    import enzomodules.problems as P
    with problems.Session(_photontest()) as s:
        g = s.grid(0)
        with P._suppress_fd_output():
            for _ in range(5):                       # heat + ionize
                s.set_boundary(0)
                s.set_dt(0, s.compute_dt(0))
                s.evolve_photons(0)
                s.advance_time(0)
            hii_before = sum(g.field("HIIDensity"))
            s.set_dt(0, s.compute_dt(0))
            energies = []
            for _ in range(5):                       # cool only
                s.solve_cooling(0)
                energies.append(sum(g.field("TotalEnergy")))
            hii_after = sum(g.field("HIIDensity"))
        assert all(b < a for a, b in zip(energies, energies[1:])), energies
        assert hii_after != hii_before               # chemistry advanced


def test_session_solve_cooling_noop():
    """Hardening: solve_cooling is a clean no-op on a problem with cooling and
    chemistry off (Toro-1), leaving the hydro state intact."""
    import enzomodules.problems as P
    with problems.Session(_toro1()) as s:
        rho0 = list(s.grid(0).field("Density"))
        with P._suppress_fd_output():
            s.set_boundary(0)
            s.set_dt(0, s.compute_dt(0))
            s.solve_cooling(0)                       # no-op, no change
        assert list(s.grid(0).field("Density")) == rho0


def test_session_init_failure_is_graceful(tmp_path):
    """Hardening: a problem whose initialization throws (ENZO_FAIL, e.g. a
    missing cooling-rate data file) raises a Python error instead of aborting
    the host process."""
    path = _param("Cooling/CoolingTest_JHW/CoolingTest_JHW.enzo")
    if not os.path.exists(path):
        pytest.skip("CoolingTest_JHW parameter file missing")
    # Run from a clean directory WITHOUT the rate-data files so init fails.
    import enzomodules.problems as P
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with pytest.raises(RuntimeError):
            with P._suppress_fd_output():
                problems.Session(path)
    finally:
        os.chdir(cwd)
    # the process survived -- a normal session still initializes fine
    with problems.Session(_toro1()) as s:
        assert s.grid(0).size > 0


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
