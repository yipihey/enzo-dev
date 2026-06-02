"""Tests for particle handling on the full grid fixture: CIC mass deposition
(legacy grid::DepositParticlePositions).

Validates the convention-independent invariants of cloud-in-cell deposition
-- the absolute mass scale follows Enzo's particle-mass unit convention, so
we check ratios and physical invariants rather than absolute values.
Requires the grid solver library (deps/build_grid.sh); skips otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)

NX = 16
H = 1.0 / NX


def _total(field, cellvol):
    return sum(field) * cellvol


def test_eight_cell_equal_split():
    """A particle equidistant from 8 cell centers deposits 1/8 to each -- the
    defining property of CIC."""
    # k*H (integer k, interior) is equidistant from 8 surrounding cell centers.
    pos = (8 * H, 8 * H, 8 * H)
    field, _ = bridge.cic_deposit([(pos[0], pos[1], pos[2], 1.0)], (NX, NX, NX))
    nz = sorted(v for v in field if abs(v) > 1e-12)
    assert len(nz) == 8
    # all eight weights equal
    assert max(nz) - min(nz) < 1e-9 * max(nz)


def test_localization():
    """One particle touches at most 8 cells (the CIC cloud)."""
    field, _ = bridge.cic_deposit([(0.515, 0.523, 0.531, 1.0)], (NX, NX, NX))
    nz = [v for v in field if abs(v) > 1e-12]
    assert 1 <= len(nz) <= 8


def test_linearity():
    """Deposited field scales linearly with particle mass."""
    f1, cv = bridge.cic_deposit([(0.515, 0.523, 0.531, 1.0)], (NX, NX, NX))
    f2, _ = bridge.cic_deposit([(0.515, 0.523, 0.531, 3.0)], (NX, NX, NX))
    for a, b in zip(f1, f2):
        if abs(a) > 1e-15:
            assert abs(b / a - 3.0) < 1e-9


def test_additivity():
    """Depositing two particles equals the sum of depositing each alone."""
    p1 = (0.30, 0.40, 0.55, 1.0)
    p2 = (0.62, 0.71, 0.48, 2.5)
    f12, cv = bridge.cic_deposit([p1, p2], (NX, NX, NX))
    f1, _ = bridge.cic_deposit([p1], (NX, NX, NX))
    f2, _ = bridge.cic_deposit([p2], (NX, NX, NX))
    for s, a, b in zip(f12, f1, f2):
        assert abs(s - (a + b)) < 1e-9 * (abs(a) + abs(b) + 1e-12)


def test_mass_conservation_independent_of_position():
    """Total deposited mass is the same wherever the particle sits."""
    totals = []
    for pos in [(0.5, 0.5, 0.5), (0.13, 0.27, 0.61), (8 * H, 8 * H, 8 * H)]:
        field, cv = bridge.cic_deposit([(pos[0], pos[1], pos[2], 1.0)], (NX, NX, NX))
        totals.append(_total(field, cv))
    assert max(totals) - min(totals) < 1e-9 * max(totals)
