"""
Sod shock tube tests.

1. Accuracy test: compare the numerical solution to the exact Riemann solution
   at a fixed resolution and check the L1 error is below a known bound.

2. Convergence rate test (Sod): the L1 error should decrease at ~O(h) as the
   resolution is refined.  This is the expected rate for a 2nd-order scheme
   applied to a problem with discontinuities (shocks and contacts reduce the
   global rate to 1st-order in L1).

3. Smooth-wave convergence test: an acoustic wave with periodic BCs and small
   amplitude returns to its initial profile after one full period.  This
   exercises the 2nd-order accuracy of MUSCL-Hancock in smooth flow: the L2
   error should decrease at O(h²).
"""

import numpy as np
import pytest

from tests.exact_riemann import exact_riemann_1d
from tests.runner import (
    BC_PERIODIC,
    BC_TRANSMISSIVE,
    DEFAULT_GAMMA,
    make_conserved,
    make_slab_grid,
    run_test,
)

# ── Sod shock-tube parameters ────────────────────────────────────────────────
SOD_RHO_L, SOD_U_L, SOD_P_L = 1.0, 0.0, 1.0
SOD_RHO_R, SOD_U_R, SOD_P_R = 0.125, 0.0, 0.1
SOD_X0 = 0.5
SOD_T = 0.2   # shock/rarefaction stays well within [0,1] at this time


def _sod_init(N, axis=0):
    """
    Sod-tube initial data on an N×4×4 slab with square cells.
    The interface is at 0.5 along `axis`.
    """
    lo, hi, res, x_cell = make_slab_grid(N, SOD_T, interface_axis=axis)
    h = 1.0 / N

    left = x_cell < SOD_X0
    rho_1d = np.where(left, SOD_RHO_L, SOD_RHO_R)
    p_1d = np.where(left, SOD_P_L, SOD_P_R)

    # Broadcast into 3-D, then place velocity along the correct axis
    def _broadcast(arr_1d):
        shape = [1, 1, 1]
        shape[axis] = N
        return arr_1d.reshape(shape) * np.ones(res)

    rho = _broadcast(rho_1d)
    p = _broadcast(p_1d)
    zeros = np.zeros(res)

    ux = _broadcast(np.where(left, SOD_U_L, SOD_U_R)) if axis == 0 else zeros
    uy = _broadcast(np.where(left, SOD_U_L, SOD_U_R)) if axis == 1 else zeros
    uz = _broadcast(np.where(left, SOD_U_L, SOD_U_R)) if axis == 2 else zeros

    data = make_conserved(rho, ux, uy, uz, p)
    return lo, hi, res, data, x_cell


def _run_sod(name, N, axis=0):
    """Run the Sod tube on an N×4×4 slab and return (x_cell, rho_profile)."""
    lo, hi, res, init_data, x_cell = _sod_init(N, axis)

    bc_lo = [BC_TRANSMISSIVE if d == axis else BC_PERIODIC for d in range(3)]
    bc_hi = [BC_TRANSMISSIVE if d == axis else BC_PERIODIC for d in range(3)]

    _, _, _, data = run_test(
        name=name,
        lo=lo, hi=hi, res=res,
        init_data=init_data,
        final_time=SOD_T,
        bc_lo=bc_lo, bc_hi=bc_hi,
    )

    # Extract 1-D density profile along the problem axis
    if axis == 0:
        rho_profile = data[:, 0, 0, 0]
    elif axis == 1:
        rho_profile = data[0, :, 0, 0]
    else:
        rho_profile = data[0, 0, :, 0]

    return x_cell, rho_profile


# ─────────────────────────────────────────────────────────────────────────────
# Test 1: single-resolution accuracy
# ─────────────────────────────────────────────────────────────────────────────

def test_sod_l1_accuracy():
    """
    At N=200, the L1 density error against the exact Riemann solution should
    be below 5e-3.  (Calibrated against the known behaviour of MUSCL-Hancock
    with HLLC on the Sod problem at this resolution and time.)
    """
    N = 200
    x_cell, rho_num = _run_sod("SodAccuracy", N)
    rho_ex, _, _ = exact_riemann_1d(
        SOD_RHO_L, SOD_U_L, SOD_P_L,
        SOD_RHO_R, SOD_U_R, SOD_P_R,
        DEFAULT_GAMMA, x_cell, SOD_T, x0=SOD_X0,
    )
    l1 = np.mean(np.abs(rho_num - rho_ex))
    print(f"\nSod accuracy (N={N}, t={SOD_T}):  L1 density error = {l1:.4e}  (threshold 5e-3)")
    assert l1 < 5e-3, f"L1 density error {l1:.4e} exceeds 5e-3 at N={N}"


# ─────────────────────────────────────────────────────────────────────────────
# Test 2: convergence rate for the Sod tube (discontinuous → ~O(h) in L1)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("N", [50, 100, 200, 400])
def test_sod_convergence(N):
    """
    Measure L1 density error at each resolution. Stored as a test attribute so
    that the caller (test_sod_convergence_rate) can aggregate them.
    """
    x_cell, rho_num = _run_sod(f"SodConvN{N}", N)
    rho_ex, _, _ = exact_riemann_1d(
        SOD_RHO_L, SOD_U_L, SOD_P_L,
        SOD_RHO_R, SOD_U_R, SOD_P_R,
        DEFAULT_GAMMA, x_cell, SOD_T, x0=SOD_X0,
    )
    l1 = np.mean(np.abs(rho_num - rho_ex))
    print(f"\nSod L1 density error:  N={N:>4}  →  {l1:.4e}")
    assert l1 > 0, "L1 error is exactly zero — something is wrong"


def test_sod_convergence_rate():
    """
    Run Sod at N=100 and N=200 and verify the L1 density error halves by at
    least a factor of 1.4 (corresponding to a convergence rate ≥ 0.5).
    The theoretical rate for a 2nd-order scheme on a problem with
    discontinuities is ~1 in L1; this test sets a conservative lower bound.
    """
    errors = {}
    for N in [100, 200]:
        x_cell, rho_num = _run_sod(f"SodRateN{N}", N)
        rho_ex, _, _ = exact_riemann_1d(
            SOD_RHO_L, SOD_U_L, SOD_P_L,
            SOD_RHO_R, SOD_U_R, SOD_P_R,
            DEFAULT_GAMMA, x_cell, SOD_T, x0=SOD_X0,
        )
        errors[N] = np.mean(np.abs(rho_num - rho_ex))

    ratio = errors[100] / errors[200]
    rate = np.log2(ratio)
    print(
        f"\nSod L1 convergence rate (100→200):  "
        f"e(100)={errors[100]:.4e}  e(200)={errors[200]:.4e}  rate={rate:.2f}  (expected ~1)"
    )
    assert rate >= 0.5, (
        f"L1 convergence rate {rate:.2f} < 0.5. "
        f"Errors: N=100 → {errors[100]:.4e}, N=200 → {errors[200]:.4e}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# Test 3: smooth-wave convergence → 2nd order in L2
# ─────────────────────────────────────────────────────────────────────────────

def _smooth_wave_init(N):
    """
    Small-amplitude right-going acoustic wave on [0,1]×slab with periodic BCs.
    After one full period t* = 1/c₀ the wave returns to its initial profile.
    """
    gamma = DEFAULT_GAMMA
    eps = 1e-3
    rho0, p0 = 1.0, 1.0
    c0 = (gamma * p0 / rho0) ** 0.5

    h = 1.0 / N
    lo = [0.0, 0.0, 0.0]
    hi = [1.0, 4 * h, 4 * h]
    res = [N, 4, 4]
    t_period = 1.0 / c0

    x_cell = np.array([lo[0] + (i + 0.5) * h for i in range(N)])
    X, _, _ = np.meshgrid(x_cell, np.zeros(4), np.zeros(4), indexing="ij")

    rho = rho0 + eps * np.sin(2 * np.pi * X)
    ux = eps * c0 * np.sin(2 * np.pi * X)
    p = p0 + eps * gamma * np.sin(2 * np.pi * X)

    data = make_conserved(rho, ux, np.zeros_like(X), np.zeros_like(X), p, gamma)
    return lo, hi, res, data, rho[:, 0, 0], t_period


@pytest.mark.parametrize("N", [32, 64, 128, 256])
def test_smooth_wave_l2_error(N):
    """Verify L2 density error is small for the acoustic wave at each resolution."""
    lo, hi, res, init_data, rho_init_1d, t_period = _smooth_wave_init(N)
    bc = [BC_PERIODIC, BC_PERIODIC, BC_PERIODIC]

    _, _, _, data = run_test(
        name=f"SmoothWaveN{N}",
        lo=lo, hi=hi, res=res,
        init_data=init_data,
        final_time=t_period,
        bc_lo=bc, bc_hi=bc,
    )

    rho_final_1d = data[:, 0, 0, 0]
    l2 = np.sqrt(np.mean((rho_final_1d - rho_init_1d) ** 2))
    print(f"\nSmooth wave L2 error vs initial condition:  N={N:>3}  →  {l2:.4e}")
    assert l2 < 0.1, f"L2 smooth-wave error {l2:.4e} suspiciously large at N={N}"


def test_smooth_wave_convergence_rate():
    """
    Self-convergence test for 2nd-order accuracy on smooth flow.

    Comparing the numerical solution to the initial condition as "exact" does
    not work here: for ε=1e-3 the accumulated nonlinear phase error (~ε² per
    period) swamps the O(h²) truncation error at all reasonable resolutions.

    Instead we use self-convergence: compare the solution on grid N to the
    block-averaged solution on grid 2N at the same final time.  The physical
    nonlinear error is identical in both solutions and cancels exactly, leaving
    only the O(h²) numerical difference.  A 2nd-order scheme gives a rate of 2.
    """
    gamma = DEFAULT_GAMMA
    c0 = (gamma * 1.0 / 1.0) ** 0.5
    t_period = 1.0 / c0
    bc = [BC_PERIODIC, BC_PERIODIC, BC_PERIODIC]

    profiles = {}
    for N in [32, 64, 128]:
        lo, hi, res, init_data, _, _ = _smooth_wave_init(N)
        _, _, _, data = run_test(
            name=f"SmoothSelfN{N}",
            lo=lo, hi=hi, res=res,
            init_data=init_data,
            final_time=t_period,
            bc_lo=bc, bc_hi=bc,
        )
        profiles[N] = data[:, 0, 0, 0]

    def block_avg(arr):
        """Average neighbouring pairs to project a fine profile onto a coarser grid."""
        return (arr[0::2] + arr[1::2]) / 2

    e1 = np.sqrt(np.mean((profiles[32] - block_avg(profiles[64])) ** 2))
    e2 = np.sqrt(np.mean((profiles[64] - block_avg(profiles[128])) ** 2))

    rate = np.log2(e1 / e2)
    print(
        f"\nSmooth wave self-convergence (L2):  "
        f"32↔64={e1:.4e}  64↔128={e2:.4e}  rate={rate:.2f}  (expected ~2)"
    )
    assert rate >= 1.6, (
        f"Self-convergence rate {rate:.2f} < 1.6 (expected ~2). "
        f"Differences: 32↔64 → {e1:.4e},  64↔128 → {e2:.4e}"
    )
