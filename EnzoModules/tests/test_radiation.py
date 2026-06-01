"""Tests that the full-grid fixture supports radiation-transport fields.

The photon-transport ray-tracer (grid::WalkPhotonPackage) is a larger
subsystem (units, HEALPix directions, species + 7 rate fields, inter-grid
transport) and is the next increment; these tests certify the foundation: a
grid carrying the radiative-transfer rate fields, and Enzo's RT field
identification (grid::IdentifyRadiativeTransferFields), work through the
fixture.  Requires the grid solver library; skips otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)


def test_rt_fields_identified_and_roundtrip():
    kph = [float(i) for i in range(20)]
    out, kphHINum, gammaNum = bridge.rt_identify(kph)
    # Enzo located both radiative-transfer rate fields.
    assert kphHINum >= 0
    assert gammaNum >= 0
    assert kphHINum != gammaNum
    # The kphHI field round-trips through the grid unchanged.
    assert out == kph


def test_rt_field_values_preserved():
    kph = [1e-12 * (i + 1) for i in range(16)]
    out, _, _ = bridge.rt_identify(kph)
    assert all(abs(a - b) <= 1e-24 for a, b in zip(out, kph))
