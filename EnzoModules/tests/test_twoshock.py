"""Replay tests for the legacy two-shock Riemann solver.

These exercise the wrapped legacy kernel against the captured golden
fixtures.  They require the compiled pilot library (build it with
``EnzoModules/deps/build_pilot.sh``); when it is absent the live-call tests
are skipped, so a checkout without the C build is still usable.
"""
import pytest

import enzomodules
from enzomodules import bridge, hydro
from enzomodules.diff import Tolerance, compare
from enzomodules.fixtures import load_dir

pytestmark = pytest.mark.skipif(
    not bridge.available(),
    reason=f"pilot library not built ({bridge.libpath()}); "
           "run EnzoModules/deps/build_pilot.sh",
)


def _fixtures():
    return load_dir(enzomodules.fixturedir("Hydro", "twoshock"))


def test_precision_contract():
    rb, ib = bridge.check_precision()
    assert rb == 8
    assert ib == 4


# Same library + same inputs -> reproduces its own captured outputs.  Allow
# only a hair of slack for the text round-trip of the fixtures.
_TOL = Tolerance(rtol=1e-13, atol=1e-300)


@pytest.mark.parametrize("fx", _fixtures(), ids=lambda f: f.name)
def test_replay_reproduces_reference(fx):
    pbar, ubar = hydro.twoshock_from_fixture(fx)
    assert compare(pbar, fx["pbar"], _TOL), f"pbar mismatch for {fx.name}"
    assert compare(ubar, fx["ubar"], _TOL), f"ubar mismatch for {fx.name}"


def test_sod_via_ergonomic_api():
    sod = {f.name: f for f in _fixtures()}["sod"]
    pbar, ubar = hydro.twoshock([1.0], [0.125], [1.0], [0.1], [0.0], [0.0],
                                gamma=1.4)
    assert compare(pbar, sod["pbar"], Tolerance(rtol=1e-13))
    assert compare(ubar, sod["ubar"], Tolerance(rtol=1e-13))
    # Physical sanity: two-shock star pressure for Sod ~ 0.3031, u* ~ 0.927.
    assert 0.29 < pbar[0] < 0.31
    assert 0.90 < ubar[0] < 0.95


def test_negative_control():
    sod = {f.name: f for f in _fixtures()}["sod"]
    pbar, _ = hydro.twoshock_from_fixture(sod)
    corrupted = list(sod["pbar"])
    corrupted[0] += 0.05
    assert not compare(pbar, corrupted, Tolerance(rtol=1e-9))
