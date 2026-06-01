"""Tests for the primordial chemistry + radiative cooling network
(legacy grid::SolveRateAndCoolEquations, the non-Grackle MultiSpecies path),
at all three MultiSpecies levels: 6 / 9 / 12 species.

Validated with conservation laws and physical direction, which hold
regardless of the exact rate values.  Requires the grid solver library
(deps/build_grid.sh); skips otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")

X = 0.76                               # hydrogen mass fraction
Y = 0.24                               # helium mass fraction
DU, LU, TU = 1.673e-24, 3.086e21, 3.156e13   # cgs: m_H, kpc, Myr
MH, KB, GAMMA, MU = 1.673e-24, 1.381e-16, 5.0 / 3.0, 1.0
TEMP_UNITS = MH * (LU / TU) ** 2 / KB


def _energy(T):
    return T / (TEMP_UNITS * (GAMMA - 1.0) * MU)


def _hydrogen_mass(sp):
    """Total H mass density across all H-bearing species (proton-mass units)."""
    h = sp["HI"][0] + sp["HII"][0]
    for k in ("HM", "H2I", "H2II"):
        if k in sp:
            h += sp[k][0]
    return h


def _initial(multispecies, T):
    """Neutral-ish gas with the H field normalized to exactly X (so H is
    conserved to X) -- cool enough that H2/HD chemistry is active."""
    sp = {
        "De": [1e-4], "HII": [X * 1e-3],
        "HeI": [Y], "HeII": [1e-10], "HeIII": [1e-10],
    }
    extra = 0.0
    if multispecies > 1:
        sp["HM"] = [1e-8]
        sp["H2I"] = [X * 1e-3]
        sp["H2II"] = [1e-10]
        extra = sp["HM"][0] + sp["H2I"][0] + sp["H2II"][0]
    # HI takes up the remainder so the H field sums to exactly X.
    sp["HI"] = [X - sp["HII"][0] - extra]
    if multispecies > 2:
        sp["DI"] = [X * 1.4e-4]
        sp["DII"] = [1e-12]
        sp["HDI"] = [1e-12]
    return sp


@pytest.mark.parametrize("multispecies,nspecies", [(1, 6), (2, 9), (3, 12)])
def test_runs_and_conserves(multispecies, nspecies):
    sp = _initial(multispecies, T=500.0)
    assert len(bridge.SPECIES_BY_LEVEL[multispecies]) == nspecies
    H0 = _hydrogen_mass(sp)
    He0 = sp["HeI"][0] + sp["HeII"][0] + sp["HeIII"][0]
    e, out = bridge.chemistry_step([1.0], [_energy(500.0)], sp, dt=0.5,
                                   multispecies=multispecies,
                                   density_units=DU, length_units=LU, time_units=TU)
    # Hydrogen and helium nuclei conserved (chemistry only changes ionization
    # / molecular state).
    assert abs(_hydrogen_mass(out) - H0) < 1e-4 * H0
    assert abs((out["HeI"][0] + out["HeII"][0] + out["HeIII"][0]) - He0) < 1e-4 * He0
    # All species stay non-negative.
    for name, vals in out.items():
        assert vals[0] >= 0.0


def test_charge_conservation_6species():
    sp = _initial(1, T=1e4)
    e, out = bridge.chemistry_step([1.0], [_energy(1e4)], sp, dt=1.0,
                                   multispecies=1,
                                   density_units=DU, length_units=LU, time_units=TU)
    charge = out["HII"][0] + out["HeII"][0] / 4.0 + out["HeIII"][0] / 2.0
    assert abs(out["De"][0] - charge) < 1e-4 * max(charge, 1e-12)


def test_h2_formation_9species():
    """At low temperature the 9-species network forms molecular hydrogen."""
    sp = _initial(2, T=300.0)
    h2_0 = sp["H2I"][0]
    e, out = bridge.chemistry_step([1.0], [_energy(300.0)], sp, dt=2.0,
                                   multispecies=2,
                                   density_units=DU, length_units=LU, time_units=TU)
    assert out["H2I"][0] > h2_0          # H2 grew


def test_hd_formation_12species():
    """The 12-species network activates deuterium chemistry (HD forms)."""
    sp = _initial(3, T=300.0)
    hd_0 = sp["HDI"][0]
    e, out = bridge.chemistry_step([1.0], [_energy(300.0)], sp, dt=2.0,
                                   multispecies=3,
                                   density_units=DU, length_units=LU, time_units=TU)
    assert out["HDI"][0] > hd_0          # HD grew from the DI reservoir
