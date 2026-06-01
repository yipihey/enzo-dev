"""Tests for the wrapped constrained-transport (CT) MHD solver,
grid::SolveMHD_Li (HydroMethod=MHD_Li, UseMHDCT=1), via grid::SolveHydroEquations
and the generic grid-fixture infrastructure.

The defining property of a CT scheme is that the discrete divergence of the
FACE-CENTERED magnetic field stays at machine precision under evolution -- the
solver advances the face B by taking the curl of an edge-centered electric
field, which is divergence-free by construction.  We set up a 1D Brio-Wu-like
MHD shock tube (gamma=2; Bx constant = 0.75, By jumping +1 -> -1 across x; the
density/pressure jump too), take a CT step, and assert:
  (a) the run succeeds,
  (b) max|divB| over active cells stays < 1e-10 both before and after,
  (c) the field actually evolves (By changes).

Requires the grid solver library, which links against the full Enzo shared
library (build with EnzoModules/deps/build_grid.sh).  Skips cleanly otherwise.
"""
import pytest

from enzomodules import bridge

pytestmark = pytest.mark.skipif(
    not bridge.grid_available(),
    reason=f"grid solver library not built ({bridge.grid_libpath()}); "
           "run EnzoModules/deps/build_grid.sh",
)

GAMMA = 2.0
NGHOST = 5          # NumberOfGhostZones the MHD_Li bridge sets
NACTIVE = 64        # active cells along x


def _brio_wu_grid():
    """Build a 1D Brio-Wu state on a 3D grid (varying along x, one active
    transverse cell).  Returns (dims, dx, fields-dict)."""
    nx = NACTIVE + 2 * NGHOST
    ny = 1 + 2 * NGHOST
    nz = 1 + 2 * NGHOST
    dims = [nx, ny, nz]
    dx = 1.0 / NACTIVE
    size = nx * ny * nz

    d = [0.0] * size
    e = [0.0] * size
    u = [0.0] * size
    v = [0.0] * size
    w = [0.0] * size
    bx = [0.0] * size
    by = [0.0] * size
    bz = [0.0] * size

    # Brio-Wu left / right states (Brio & Wu 1988).
    Bx = 0.75
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = i + nx * (j + ny * k)
                # discontinuity at the active midpoint
                left = (i - NGHOST) < (NACTIVE // 2)
                rho = 1.0 if left else 0.125
                p = 1.0 if left else 0.1
                By = 1.0 if left else -1.0
                d[idx] = rho
                u[idx] = 0.0
                v[idx] = 0.0
                w[idx] = 0.0
                bx[idx] = Bx
                by[idx] = By
                bz[idx] = 0.0
                # specific total energy = internal + 0.5 v^2 + 0.5 B^2 / rho
                eint = p / ((GAMMA - 1.0) * rho)
                emag = 0.5 * (Bx * Bx + By * By) / rho
                e[idx] = eint + emag
    return dims, dx, dict(d=d, e=e, u=u, v=v, w=w, bx=bx, by=by, bz=bz)


def test_mhdct_step_runs_and_preserves_divb():
    """One CT step on a Brio-Wu tube: succeeds, keeps div B ~ 0, evolves By."""
    dims, dx, f = _brio_wu_grid()
    dt = 0.2 * dx  # CFL-safe for fast speed ~ a few here

    by0 = list(f["by"])
    (d, e, u, v, w, bx, by, bz,
     divb_before, divb_after) = bridge.mhdct_step(
        dims, f["d"], f["e"], f["u"], f["v"], f["w"],
        f["bx"], f["by"], f["bz"], dx, dt, GAMMA)

    # (b) the CT guarantee: div B at machine precision both before & after.
    assert divb_before < 1e-10, f"initial max|divB| = {divb_before:.3e}"
    assert divb_after < 1e-10, f"post-step max|divB| = {divb_after:.3e}"

    # (c) the field actually evolved.
    dby = max(abs(a - b) for a, b in zip(by, by0))
    assert dby > 1e-6, f"By did not change (max |dBy| = {dby:.3e})"

    # sanity: density stays positive and finite.
    assert all(x == x and x > 0 for x in d)


def test_mhdct_multistep_divb_stays_zero():
    """A few sequential CT steps keep div B at machine precision (the invariant
    holds under repeated evolution, not just one step).

    The steps are taken on a single persistent grid (nsteps=4) so the
    face-centered field -- the true CT state -- carries over natively; the
    reported max|divB| is the worst over all four steps.
    """
    dims, dx, f = _brio_wu_grid()
    dt = 0.2 * dx
    by0 = list(f["by"])

    (d, e, u, v, w, bx, by, bz,
     divb_before, divb_after) = bridge.mhdct_step(
        dims, f["d"], f["e"], f["u"], f["v"], f["w"],
        f["bx"], f["by"], f["bz"], dx, dt, GAMMA, nsteps=4)

    assert divb_before < 1e-10, f"initial max|divB| = {divb_before:.3e}"
    assert divb_after < 1e-10, f"max|divB| over 4 CT steps = {divb_after:.3e}"
    dby = max(abs(a - b) for a, b in zip(by, by0))
    assert dby > 1e-6, f"By did not evolve over 4 steps (max |dBy| = {dby:.3e})"
