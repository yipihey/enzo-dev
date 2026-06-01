"""Tests for the AMR control path: cell flagging (grid::SetFlaggingField) and
clustering flagged cells into child grids (Berger-Rigoutsos via ProtoSubgrid +
IdentifyNewSubgridsBySignature).

Requires the grid solver library (deps/build_grid.sh); skips otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")

NG = 3


def test_flag_density_jump():
    """Slope flagging marks the cells straddling a sharp density jump and
    leaves the flat regions alone."""
    nx = 32
    dimx = nx + 2 * NG
    jump = nx // 2
    density = [1.0 if (i - NG) < jump else 0.1 for i in range(dimx)]
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), density,
                                        slope_threshold=0.3)
    flagged = [i - NG for i in range(NG, dimx - NG) if flagging[i]]
    assert count > 0
    # flagged cells are exactly at the discontinuity...
    assert all(abs(c - jump) <= 1 for c in flagged)
    # ...and the far-field flat regions are not flagged.
    assert flagging[NG + 2] == 0 and flagging[dimx - NG - 3] == 0


def test_flag_uniform_nothing():
    """A uniform field flags no cells."""
    nx = 16
    dimx = nx + 2 * NG
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), [1.0] * dimx,
                                        slope_threshold=0.3)
    assert count == 0
    assert sum(flagging) == 0


def _block(fl, dimx, i0, i1, j0, j1):
    for j in range(j0, j1):
        for i in range(i0, i1):
            fl[(i + NG) + (j + NG) * dimx] = 1


def test_cluster_single_block():
    """One flagged block clusters into one subgrid tightly enclosing it."""
    nx = ny = 16
    dimx = dimy = nx + 2 * NG
    fl = [0] * (dimx * dimy)
    _block(fl, dimx, 4, 10, 5, 9)
    boxes = bridge.cluster(2, (dimx, dimy, 1), fl)
    assert len(boxes) == 1
    (xlo, xhi), (ylo, yhi), _ = boxes[0]
    assert (xlo, xhi) == (4, 10)
    assert (ylo, yhi) == (5, 9)


def test_cluster_two_blocks():
    """Two separated blocks cluster into two disjoint subgrids."""
    nx = ny = 16
    dimx = dimy = nx + 2 * NG
    fl = [0] * (dimx * dimy)
    _block(fl, dimx, 2, 5, 2, 5)
    _block(fl, dimx, 10, 14, 10, 14)
    boxes = bridge.cluster(2, (dimx, dimy, 1), fl)
    assert len(boxes) == 2
    # each block is enclosed by exactly one subgrid; boxes are disjoint in x
    xranges = sorted((b[0][0], b[0][1]) for b in boxes)
    assert xranges[0] == (2, 5)
    assert xranges[1] == (10, 14)


def test_cluster_efficiency():
    """The clustered subgrid is minimal -- it contains no slack beyond the
    flagged region (every boundary row/col is flagged)."""
    nx = ny = 24
    dimx = dimy = nx + 2 * NG
    fl = [0] * (dimx * dimy)
    _block(fl, dimx, 6, 18, 8, 16)
    boxes = bridge.cluster(2, (dimx, dimy, 1), fl)
    assert len(boxes) == 1
    (xlo, xhi), (ylo, yhi), _ = boxes[0]
    assert (xlo, xhi, ylo, yhi) == (6, 18, 8, 16)
