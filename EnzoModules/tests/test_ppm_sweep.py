"""Tests for the composite PPM 1D sweep (inteuler -> twoshock ->
flux_twoshock -> euler).

Replay tests reproduce captured single-step outputs bit-for-bit; the
integration test evolves a Sod tube and checks it against the exact Riemann
solution.  All require the compiled pilot library and skip without it.
"""
import pytest

import enzomodules
from enzomodules import bridge, hydro
from enzomodules.diff import Tolerance, compare
from enzomodules.examples import ppm_sod, riemann
from enzomodules.fixtures import load_dir

pytestmark = pytest.mark.skipif(
    not bridge.available(),
    reason=f"pilot library not built ({bridge.libpath()}); "
           "run EnzoModules/deps/build_pilot.sh",
)

_STATE_FIELDS = ["dslice_out", "eslice_out", "uslice_out", "vslice_out", "wslice_out"]


def _fixtures():
    return load_dir(enzomodules.fixturedir("Hydro", "ppm_sweep_1d"))


# Same library + same inputs -> reproduces its own captured outputs.
_TOL = Tolerance(rtol=1e-12, atol=1e-290)


@pytest.mark.parametrize("fx", _fixtures(), ids=lambda f: f.name)
def test_single_step_replay(fx):
    d, e, u, v, w = hydro.ppm_sweep_from_fixture(fx)
    got = {"dslice_out": d, "eslice_out": e, "uslice_out": u,
           "vslice_out": v, "wslice_out": w}
    for field in _STATE_FIELDS:
        assert compare(got[field], fx[field], _TOL), f"{field} mismatch in {fx.name}"


def test_negative_control():
    fx = _fixtures()[0]
    d, *_ = hydro.ppm_sweep_from_fixture(fx)
    corrupted = list(fx["dslice_out"])
    corrupted[len(corrupted) // 2] += 0.1
    assert not compare(d, corrupted, Tolerance(rtol=1e-9))


def test_sod_matches_exact_riemann():
    """End-to-end: evolve Sod with the wrapped PPM solver, compare to truth."""
    gamma = 1.4
    g, nsteps = ppm_sod.run(t_final=0.2, nx=200, gamma=gamma)
    x = g.x_centers()
    rho = g.active(g.d)
    u = g.active(g.u)
    p = g.active(g.pressure())

    rho_ex, u_ex, p_ex = riemann.sample((1.0, 0.0, 1.0), (0.125, 0.0, 0.1),
                                        gamma, 0.5, x, 0.2)
    n = len(rho)
    l1_rho = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / n
    l1_p = sum(abs(a - b) for a, b in zip(p, p_ex)) / n

    assert nsteps > 50                    # actually integrated in time
    assert l1_rho < 1.5e-2, f"L1 density error {l1_rho:.3e}"
    assert l1_p < 1.5e-2, f"L1 pressure error {l1_p:.3e}"
    # Plateaus pinned by the analytic solution.
    assert abs(max(rho) - 1.0) < 1e-3
    assert abs(min(rho) - 0.125) < 1e-3


def test_star_state_sanity():
    """The exact solver used for validation reproduces the classic Sod star
    pressure/velocity (guards the validator itself)."""
    pstar, ustar = riemann.star_state((1.0, 0.0, 1.0), (0.125, 0.0, 0.1), 1.4)
    assert abs(pstar - 0.30313) < 1e-3
    assert abs(ustar - 0.92745) < 1e-3
