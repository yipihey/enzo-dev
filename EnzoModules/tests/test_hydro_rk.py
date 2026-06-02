"""Tests for the wrapped hydro_rk (RK MUSCL) solvers, hydro + Dedner MHD.

These require the hydro_rk shared library, which links against the full Enzo
shared library (build with EnzoModules/deps/build_hydro_rk.sh).  They skip
cleanly when it is absent.
"""
import pytest

from enzomodules import bridge
from enzomodules.examples import hydro_rk_sod, mhd_brio_wu, riemann

pytestmark = pytest.mark.skipif(
    not bridge.hydrork_available(),
    reason=f"hydro_rk library not built ({bridge.hydrork_libpath()}); "
           "run EnzoModules/deps/build_hydro_rk.sh",
)


@pytest.mark.parametrize("solver", ["hll", "hllc", "llf"])
def test_hydro_sod_matches_exact(solver):
    """Each hydro Riemann solver evolves Sod to match the exact solution."""
    x, rho, u, p, nsteps = hydro_rk_sod.run(solver=solver, t_final=0.2, nx=200)
    rho_ex, u_ex, p_ex = riemann.sample((1.0, 0.0, 1.0), (0.125, 0.0, 0.1),
                                        1.4, 0.5, x, 0.2)
    n = len(rho)
    l1 = sum(abs(a - b) for a, b in zip(rho, rho_ex)) / n
    assert nsteps > 50
    assert l1 < 6e-3, f"{solver} L1 density error {l1:.3e}"
    assert abs(max(rho) - 1.0) < 1e-3
    assert abs(min(rho) - 0.125) < 1e-3


def test_hydro_flux_sanity():
    """At a t=0 Sod interface, mass flux is zero away from the jump and the
    momentum flux equals the (static) pressure."""
    nghost = 3
    ncells = 10 + 2 * nghost
    gamma = 1.4
    rows = [[0.0] * ncells for _ in range(5)]
    for c in range(ncells):
        left = (c - nghost) < 5
        rho, pr = (1.0, 1.0) if left else (0.125, 0.1)
        rows[0][c] = rho
        rows[1][c] = pr / ((gamma - 1.0) * rho)   # internal energy
    flux = bridge.hydro_rk_line(bridge.HLL, rows, gamma=gamma)
    # momentum-x flux equals pressure in the uniform regions
    assert abs(flux[1][0] - 1.0) < 1e-12
    assert abs(flux[1][-1] - 0.1) < 1e-12
    # mass flux vanishes in the uniform regions, non-zero only at the jump
    assert abs(flux[0][0]) < 1e-12 and abs(flux[0][-1]) < 1e-12
    assert any(abs(flux[0][i]) > 1e-3 for i in range(len(flux[0])))


def test_dedner_mhd_brio_wu():
    """The Dedner HLLD solver reproduces the Brio-Wu shock tube and keeps the
    normal field Bx divergence-free (constant)."""
    x, rho, vx, By, Bx, nsteps = mhd_brio_wu.run(solver="hlld", t_final=0.1,
                                                 nx=400)
    assert nsteps > 50
    # Known Brio-Wu envelope (gamma=2).
    assert 0.10 < min(rho) < 0.16
    assert abs(max(rho) - 1.0) < 1e-2
    assert min(By) < -0.5 and max(By) > 0.5
    # Dedner cleaning: the normal field stays at its initial constant value.
    assert all(abs(b - mhd_brio_wu.BRIO_WU_BX) < 1e-6 for b in Bx)


def test_dedner_coupling_present():
    """The GLM cleaning couples Bx and Phi through C_h: a jump in the normal
    field Bx (a divergence error) produces a Bx flux = -0.5*C_h*dBx, so it
    scales linearly with C_h."""
    nghost = 3
    ncells = 8 + 2 * nghost
    rows = [[0.0] * ncells for _ in range(9)]
    for c in range(ncells):
        rows[0][c] = 1.0           # rho
        rows[1][c] = 1.0           # eint
        rows[5][c] = 1.0 if (c - nghost) < 4 else -1.0   # Bx jump (divergence)
        # Phi = 0 everywhere
    f_lo = bridge.mhd_rk_line(bridge.HLLD, rows, c_h=1.0, gamma=2.0)
    f_hi = bridge.mhd_rk_line(bridge.HLLD, rows, c_h=2.0, gamma=2.0)
    jump = len(f_lo[5]) // 2
    assert abs(f_lo[5][jump]) > 1e-6
    # Doubling C_h doubles the cleaning flux (within reconstruction effects).
    assert abs(f_hi[5][jump]) > 1.5 * abs(f_lo[5][jump])
