"""Tests for the AMR hierarchy inter-grid operators that couple refinement
levels, certified by conservation / accuracy:

  1. RESTRICTION  -- grid::ProjectSolutionToParentGrid (fine -> coarse average)
  2. PROLONGATION -- grid::InterpolateFieldValues      (coarse -> fine interp)
  3. REFLUXING    -- grid::CorrectForRefinedFluxes     (flux correction)

A parent grid (N active cells over the unit domain) and a child grid covering
parent active cells [N/4, 3N/4) at 2x resolution (refinement factor 2).

Requires the grid solver library (deps/build_grid.sh); skips otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh")

NG = 3


# ---- 1. RESTRICTION: conservation (volume average) -----------------------

def test_restriction_linear_ramp_is_volume_average():
    """Fill the child with a linear ramp f = a + b*x at child cell centers,
    project down to the parent, and assert each overlapped parent cell equals
    the mean of the 2 child cells it covers (the conservation property)."""
    n = 16
    a, b = 0.7, 2.3
    dxp = 1.0 / n
    cdx = dxp / 2.0
    cleft = (n // 4) * dxp
    cdim = n + 2 * NG

    # child cell centers (active cell j has center cleft + (j+0.5)*cdx);
    # ghost cells continue the same ramp so nothing is special at the edges.
    child = []
    for i in range(cdim):
        j = i - NG                       # active-cell index (negative in ghosts)
        x = cleft + (j + 0.5) * cdx
        child.append(a + b * x)

    parent, pstart, noverlap, ng = bridge.project_to_parent(n, child)
    assert ng == NG
    assert noverlap == n // 2

    max_err = 0.0
    for p in range(noverlap):
        cl = 2 * p              # child active indices covered by parent cell
        cr = 2 * p + 1
        expected = 0.5 * (child[NG + cl] + child[NG + cr])
        got = parent[NG + pstart + p]
        max_err = max(max_err, abs(got - expected))
    assert max_err < 1e-12, f"restriction conservation error {max_err:.3e}"


def test_restriction_uniform_field_preserved():
    """A uniform child field projects to the same uniform parent value."""
    n = 16
    cdim = n + 2 * NG
    child = [3.14159] * cdim
    parent, pstart, noverlap, ng = bridge.project_to_parent(n, child)
    for p in range(noverlap):
        assert abs(parent[NG + pstart + p] - 3.14159) < 1e-12


# ---- 2. PROLONGATION: accuracy (linear exact) ----------------------------

def test_prolongation_linear_is_exact():
    """Fill the parent with a linear field (exact for SecondOrderA), interpolate
    up to the child, and assert the child cell-center values match the analytic
    field to ~1e-10 (interpolation accuracy)."""
    n = 16
    a, b = 1.5, 0.8
    dxp = 1.0 / n
    pdim = n + 2 * NG

    # parent cell centers; ghosts continue the ramp so the stencil is clean.
    parent = []
    for i in range(pdim):
        j = i - NG
        x = (j + 0.5) * dxp
        parent.append(a + b * x)

    child, cactive, ng, cleft, cdx = bridge.interpolate_to_child(n, parent)
    assert ng == NG and cactive == n

    max_err = 0.0
    for j in range(cactive):
        xc = cleft + (j + 0.5) * cdx
        expected = a + b * xc
        got = child[NG + j]
        max_err = max(max_err, abs(got - expected))
    assert max_err < 1e-10, f"prolongation accuracy error {max_err:.3e}"


def test_prolongation_conservation_child_average_equals_parent():
    """Conservative interpolation: the mean of the 2 child cells under each
    parent cell equals that parent cell's value (for a linear field)."""
    n = 16
    a, b = 2.0, -0.5
    dxp = 1.0 / n
    pdim = n + 2 * NG
    parent = [a + b * ((i - NG) + 0.5) * dxp for i in range(pdim)]

    child, cactive, ng, cleft, cdx = bridge.interpolate_to_child(n, parent)
    pstart = n // 4
    max_err = 0.0
    for p in range(cactive // 2):
        avg = 0.5 * (child[NG + 2 * p] + child[NG + 2 * p + 1])
        pval = parent[NG + pstart + p]
        max_err = max(max_err, abs(avg - pval))
    assert max_err < 1e-10, f"prolongation conservation error {max_err:.3e}"


# ---- 3. REFLUXING: conservative flux correction --------------------------

def test_refluxing_identity_noop():
    """Equal initial and refined fluxes leave the field unchanged (no-op)."""
    n = 16
    gdim = n + 2 * NG
    density = [1.0 + 0.01 * i for i in range(gdim)]
    lo, hi = n // 4, 3 * n // 4 - 1
    out, lc, rc = bridge.correct_refined_fluxes(
        n, density, lo, hi,
        init_left=0.42, refined_left=0.42,
        init_right=-0.17, refined_right=-0.17)
    max_err = max(abs(out[i] - density[i]) for i in range(gdim))
    assert max_err < 1e-12, f"reflux identity drifted {max_err:.3e}"


def test_refluxing_known_difference():
    """A known flux difference corrects exactly the two coarse cells adjacent
    to the fine boundary by (Initial - Refined):
        left  cell:  d += (init_left  - refined_left)
        right cell:  d -= (init_right - refined_right)
    and leaves all other cells untouched."""
    n = 16
    gdim = n + 2 * NG
    density = [5.0] * gdim
    lo, hi = n // 4, 3 * n // 4 - 1

    il, rl = 0.30, 0.10            # left face: init, refined
    ir, rr = 0.50, 0.05            # right face: init, refined
    out, lc, rcell = bridge.correct_refined_fluxes(
        n, density, lo, hi, il, rl, ir, rr)

    # corrected coarse cells are the active cells just outside the subgrid.
    assert lc == lo - 1 and rcell == hi + 1
    left_idx = NG + lc
    right_idx = NG + rcell
    exp_left = 5.0 + (il - rl)
    exp_right = 5.0 - (ir - rr)
    assert abs(out[left_idx] - exp_left) < 1e-6, \
        f"left correction {out[left_idx]} != {exp_left}"
    assert abs(out[right_idx] - exp_right) < 1e-6, \
        f"right correction {out[right_idx]} != {exp_right}"

    # every other cell unchanged.
    for i in range(gdim):
        if i in (left_idx, right_idx):
            continue
        assert abs(out[i] - 5.0) < 1e-6, f"cell {i} drifted to {out[i]}"
