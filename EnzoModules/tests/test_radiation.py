"""Tests for radiative transfer: field plumbing and the photon-transport
ray-tracer.

- Field support: a grid carrying the radiative-transfer rate fields and Enzo's
  field identification (grid::IdentifyRadiativeTransferFields).
- Ray tracing: a single photon package fired through a uniform-HI grid
  (grid::WalkPhotonPackage) must obey Beer-Lambert attenuation,
  N(L) = N0 * exp(-n_HI * sigma * L), with sigma from Enzo's own
  FindCrossSection.

Requires the grid solver library; skips otherwise.
"""
import math

import pytest

from enzomodules import bridge

MH = 1.673e-24
DU, LU, TU = 1.673e-24, 3.086e21, 3.156e13   # m_H, kpc, Myr

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


def test_hi_cross_section_threshold():
    """The HI cross-section just above threshold is ~6e-18 cm^2 and falls with
    energy."""
    sig14 = bridge.hi_cross_section(14.0)
    sig50 = bridge.hi_cross_section(50.0)
    assert 5e-18 < sig14 < 7e-18
    assert sig50 < sig14            # photo-ionization cross-section decreases


@pytest.mark.parametrize("n_HI", [3e-5, 1e-4, 3e-4, 1e-3])
def test_beer_lambert(n_HI):
    """A photon ray through uniform HI attenuates as exp(-n*sigma*L)."""
    energy = 14.0
    N0 = 1e50
    sigma = bridge.hi_cross_section(energy)
    photons_final, radius, kph_sum = bridge.raytrace_uniform(
        n_HI, energy=energy, photons=N0, path_fraction=0.3,
        density_units=DU, length_units=LU, time_units=TU)
    assert photons_final > 0          # not fully absorbed at these depths
    tau = (n_HI * DU / MH) * sigma * (radius * LU)
    expected = N0 * math.exp(-tau)
    assert abs(photons_final - expected) <= 5e-3 * expected, \
        f"N/N0={photons_final/N0:.4e} vs exp(-tau)={expected/N0:.4e} (tau={tau:.3f})"


def test_optically_thick_absorbs_everything():
    """A very optically-thick ray is fully absorbed (photon deleted)."""
    photons_final, radius, kph_sum = bridge.raytrace_uniform(
        1e-2, energy=14.0, photons=1e50, path_fraction=0.3,
        density_units=DU, length_units=LU, time_units=TU)
    # Either driven to ~0 or flagged deleted (negative sentinel).
    assert photons_final <= 1e50 * 1e-4
    assert kph_sum > 0                # photons were deposited as ionizations
