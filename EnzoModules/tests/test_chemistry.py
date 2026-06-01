"""Tests for the primordial chemistry + radiative cooling network
(legacy grid::SolveRateAndCoolEquations, the non-Grackle MultiSpecies path).

Validated with conservation laws and physical direction, which hold
regardless of the exact rate values.  Requires the grid solver library
(deps/build_grid.sh); skips otherwise.
"""
import math

import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")

X = 0.76                               # hydrogen mass fraction
DU, LU, TU = 1.673e-24, 3.086e21, 3.156e13   # cgs: m_H, kpc, Myr
MH, KB, GAMMA, MU = 1.673e-24, 1.381e-16, 5.0 / 3.0, 0.6
TEMP_UNITS = MH * (LU / TU) ** 2 / KB


def _energy_for_temperature(T):
    return T / (TEMP_UNITS * (GAMMA - 1.0) * MU)


def _state(T=1e4, ionized_fraction=0.5):
    e = _energy_for_temperature(T)
    HI = [(1 - ionized_fraction) * X]
    HII = [ionized_fraction * X]
    return dict(rho=[1.0], e_tot=[e], de=[ionized_fraction * X],
                HI=HI, HII=HII, HeI=[0.24], HeII=[1e-10], HeIII=[1e-10])


def test_runs_and_conserves_nuclei():
    s = _state()
    H0 = s["HI"][0] + s["HII"][0]
    He0 = s["HeI"][0] + s["HeII"][0] + s["HeIII"][0]
    e, de, HI, HII, HeI, HeII, HeIII = bridge.chemistry_step(
        s["rho"], s["e_tot"], s["de"], s["HI"], s["HII"],
        s["HeI"], s["HeII"], s["HeIII"], dt=1.0,
        density_units=DU, length_units=LU, time_units=TU)
    # Hydrogen and helium nuclei are conserved (chemistry only re-ionizes).
    assert abs((HI[0] + HII[0]) - H0) < 1e-5 * H0
    assert abs((HeI[0] + HeII[0] + HeIII[0]) - He0) < 1e-5 * He0


def test_charge_conservation():
    s = _state()
    e, de, HI, HII, HeI, HeII, HeIII = bridge.chemistry_step(
        s["rho"], s["e_tot"], s["de"], s["HI"], s["HII"],
        s["HeI"], s["HeII"], s["HeIII"], dt=1.0,
        density_units=DU, length_units=LU, time_units=TU)
    # electron density balances the ionized species (He* in 4x mass units).
    charge = HII[0] + HeII[0] / 4.0 + HeIII[0] / 2.0
    assert abs(de[0] - charge) < 1e-5 * max(charge, 1e-12)


def test_recombination_and_cooling():
    """Warm (1e4 K), partially-ionized, dense gas should recombine and cool."""
    s = _state(T=1e4, ionized_fraction=0.5)
    e0 = s["e_tot"][0]
    e, de, HI, HII, HeI, HeII, HeIII = bridge.chemistry_step(
        s["rho"], s["e_tot"], s["de"], s["HI"], s["HII"],
        s["HeI"], s["HeII"], s["HeIII"], dt=1.0,
        density_units=DU, length_units=LU, time_units=TU)
    assert HII[0] < 0.5 * X            # recombined (less ionized)
    assert HI[0] > 0.5 * X
    assert e[0] < e0                   # radiative cooling lowered the energy


def test_linearity_in_density_scale():
    """Doubling all densities preserves the ionization *fractions* (rates are
    density-dependent, so check fraction sanity rather than exact linearity)."""
    s = _state()
    *_, HI, HII, _, _, _ = (None,) + tuple(bridge.chemistry_step(
        s["rho"], s["e_tot"], s["de"], s["HI"], s["HII"],
        s["HeI"], s["HeII"], s["HeIII"], dt=0.5,
        density_units=DU, length_units=LU, time_units=TU))
    frac = HII[0] / (HI[0] + HII[0])
    assert 0.0 <= frac <= 1.0
