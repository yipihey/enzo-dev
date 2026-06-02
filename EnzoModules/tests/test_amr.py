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


def test_flag_density_jump_slope():
    """Slope flagging (method 1) marks the cells straddling a sharp density
    jump and leaves the flat regions alone."""
    nx = 32
    dimx = nx + 2 * NG
    jump = nx // 2
    density = [1.0 if (i - NG) < jump else 0.1 for i in range(dimx)]
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), density,
                                        method=bridge.FLAG_SLOPE, threshold=0.3)
    flagged = [i - NG for i in range(NG, dimx - NG) if flagging[i]]
    assert count > 0
    assert all(abs(c - jump) <= 1 for c in flagged)
    assert flagging[NG + 2] == 0 and flagging[dimx - NG - 3] == 0


def test_flag_uniform_nothing():
    """A uniform field flags no cells (slope)."""
    nx = 16
    dimx = nx + 2 * NG
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), [1.0] * dimx,
                                        method=bridge.FLAG_SLOPE, threshold=0.3)
    assert count == 0
    assert sum(flagging) == 0


def test_flag_overdensity_mass():
    """Mass/overdensity flagging (method 2) marks the cells in a dense blob."""
    nx = 32
    dimx = nx + 2 * NG
    dx = 1.0 / nx
    density = [5.0 if abs((i - NG) - 16) < 3 else 1.0 for i in range(dimx)]
    # cell mass = density*dx; flag cells with density > 3 -> threshold 3*dx.
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), density,
                                        method=bridge.FLAG_MASS,
                                        threshold=3.0 * dx, dx=dx)
    flagged = [i - NG for i in range(NG, dimx - NG) if flagging[i]]
    assert flagged == [14, 15, 16, 17, 18]


def test_flag_second_derivative():
    """Second-derivative flagging (method 15) marks the high-curvature region
    of a smooth bump (and not the flat far field)."""
    import math
    nx = 32
    dimx = nx + 2 * NG
    dx = 1.0 / nx
    density = [1.0 + math.exp(-((i - NG - 16) / 2.0) ** 2) for i in range(dimx)]
    flagging, count = bridge.flag_cells(1, (dimx, 1, 1), density,
                                        method=bridge.FLAG_SECOND_DERIVATIVE,
                                        threshold=0.05, dx=dx)
    flagged = [i - NG for i in range(NG, dimx - NG) if flagging[i]]
    assert count > 0
    assert all(abs(c - 16) <= 8 for c in flagged)     # localized at the bump
    assert flagging[NG] == 0 and flagging[dimx - NG - 1] == 0


def test_flag_shear():
    """Shear flagging (method 9) marks the velocity shear interface."""
    nx = 32
    ny = 8
    dimx = nx + 2 * NG
    dimy = ny + 2 * NG
    dx = 1.0 / nx
    size = dimx * dimy
    density = [1.0] * size
    vy = [0.0] * size
    for j in range(dimy):
        for i in range(dimx):
            vy[i + j * dimx] = 1.0 if (i - NG) < nx // 2 else -1.0
    flagging, count = bridge.flag_cells(2, (dimx, dimy, 1), density,
                                        method=bridge.FLAG_SHEAR,
                                        threshold=0.5, v=vy, dx=dx)
    cols = sorted(set(i - NG for j in range(NG, dimy - NG)
                      for i in range(NG, dimx - NG) if flagging[i + j * dimx]))
    assert count > 0
    assert all(abs(c - nx // 2) <= 1 for c in cols)   # at the shear interface


def test_flag_must_refine_region():
    """Geometric must-refine-region flagging (method 12, 3D) flags exactly the
    cells inside the requested box."""
    n = 12
    dim = n + 2 * NG
    flagging, count = bridge.flag_region((dim, dim, dim),
                                         (0.25, 0.25, 0.25), (0.5, 0.5, 0.5),
                                         dx=1.0 / n)
    # box [0.25,0.5] at n=12 -> active cells 3,4,5 in each dim -> 3^3 = 27
    assert count == 27
    xs = sorted(set(i - NG for k in range(NG, dim - NG)
                    for j in range(NG, dim - NG)
                    for i in range(NG, dim - NG)
                    if flagging[i + dim * (j + dim * k)]))
    assert xs == [3, 4, 5]


def test_flag_jeans():
    """Jeans-length flagging (method 6): a dense clump is under-resolved and
    flagged; diffuse gas is well-resolved and not flagged."""
    import math
    nx = 32
    dimx = nx + 2 * NG
    DU, LU, TU = 1.673e-24, 3.086e21, 3.156e13
    mh, kB, gamma, mu = 1.673e-24, 1.381e-16, 5.0 / 3.0, 0.6
    temp_units = mh * (LU / TU) ** 2 / kB
    eint = 50.0 / (temp_units * (gamma - 1) * mu)
    energy = [eint] * dimx
    clump = [1e4 if abs((i - NG) - 16) < 4 else 1.0 for i in range(dimx)]
    fl, count = bridge.flag_jeans(1, (dimx, 1, 1), clump, energy,
                                  density_units=DU, length_units=LU,
                                  time_units=TU, dx=1.0 / nx)
    flagged = [i - NG for i in range(NG, dimx - NG) if fl[i]]
    assert count > 0
    assert all(13 <= c <= 19 for c in flagged)        # the dense clump
    # diffuse uniform gas: nothing under-resolved
    fl2, count2 = bridge.flag_jeans(1, (dimx, 1, 1), [1.0] * dimx, energy,
                                    density_units=DU, length_units=LU,
                                    time_units=TU, dx=1.0 / nx)
    assert count2 == 0


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
