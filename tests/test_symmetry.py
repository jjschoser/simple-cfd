"""
Symmetry and rotational-invariance tests.

The Euler equations are rotationally symmetric, and a well-implemented
Cartesian finite volume scheme should preserve this symmetry exactly (up to
floating-point rounding) when the grid and initial conditions are related by
rotation.

We verify this by running the same Sod shock tube in the x, y, and z
directions on identical thin-slab grids (with square cells) and comparing
the resulting 1-D density profiles.  The profiles must agree to within a
tight relative tolerance.
"""

import numpy as np
import pytest

from tests.runner import (
    BC_PERIODIC,
    BC_TRANSMISSIVE,
    DEFAULT_GAMMA,
    make_conserved,
    make_slab_grid,
    run_test,
)

SOD_RHO_L, SOD_U_L, SOD_P_L = 1.0, 0.0, 1.0
SOD_RHO_R, SOD_U_R, SOD_P_R = 0.125, 0.0, 0.1
SOD_X0 = 0.5
SOD_T = 0.2
N_SYM = 100


def _run_sod_axis(axis, suffix=""):
    """
    Run the Sod tube along `axis` on an N_SYM×4×4 slab.
    Returns the 1-D density profile along that axis.
    """
    lo, hi, res, x_cell = make_slab_grid(N_SYM, SOD_T, interface_axis=axis)
    h = 1.0 / N_SYM

    left = x_cell < SOD_X0
    rho_1d = np.where(left, SOD_RHO_L, SOD_RHO_R)
    p_1d = np.where(left, SOD_P_L, SOD_P_R)

    def _bc(arr_1d):
        shape = [1, 1, 1]
        shape[axis] = N_SYM
        return arr_1d.reshape(shape) * np.ones(res)

    rho = _bc(rho_1d)
    p = _bc(p_1d)
    zeros = np.zeros(res)

    # Velocity along the shock-tube axis only
    u_1d = np.where(left, SOD_U_L, SOD_U_R)
    ux = _bc(u_1d) if axis == 0 else zeros
    uy = _bc(u_1d) if axis == 1 else zeros
    uz = _bc(u_1d) if axis == 2 else zeros

    init_data = make_conserved(rho, ux, uy, uz, p)

    bc_lo = [BC_TRANSMISSIVE if d == axis else BC_PERIODIC for d in range(3)]
    bc_hi = [BC_TRANSMISSIVE if d == axis else BC_PERIODIC for d in range(3)]

    axis_label = ["X", "Y", "Z"][axis]
    _, _, _, data = run_test(
        name=f"SymSod{axis_label}{suffix}",
        lo=lo, hi=hi, res=res,
        init_data=init_data,
        final_time=SOD_T,
        bc_lo=bc_lo,
        bc_hi=bc_hi,
    )

    if axis == 0:
        return data[:, 0, 0, 0]
    elif axis == 1:
        return data[0, :, 0, 0]
    else:
        return data[0, 0, :, 0]


def test_symmetry_x_vs_y():
    """Sod tube along x and along y must produce the same density profile."""
    rho_x = _run_sod_axis(0, suffix="vs_y_x")
    rho_y = _run_sod_axis(1, suffix="vs_y_y")

    max_diff = np.max(np.abs(rho_x - rho_y))
    mean_rho = 0.5 * (np.mean(rho_x) + np.mean(rho_y))
    rel_diff = max_diff / mean_rho

    print(f"\nSymmetry x↔y:  max |Δρ| = {max_diff:.2e}  relative = {rel_diff:.2e}  (threshold 1e-10)")
    assert rel_diff < 1e-10, (
        f"x- vs y-directed Sod profiles differ by {rel_diff:.2e} (max abs: {max_diff:.2e}). "
        "The scheme is not rotationally symmetric in x-y."
    )


def test_symmetry_x_vs_z():
    """Sod tube along x and along z must produce the same density profile."""
    rho_x = _run_sod_axis(0, suffix="vs_z_x")
    rho_z = _run_sod_axis(2, suffix="vs_z_z")

    max_diff = np.max(np.abs(rho_x - rho_z))
    mean_rho = 0.5 * (np.mean(rho_x) + np.mean(rho_z))
    rel_diff = max_diff / mean_rho

    print(f"\nSymmetry x↔z:  max |Δρ| = {max_diff:.2e}  relative = {rel_diff:.2e}  (threshold 1e-10)")
    assert rel_diff < 1e-10, (
        f"x- vs z-directed Sod profiles differ by {rel_diff:.2e} (max abs: {max_diff:.2e}). "
        "The scheme is not rotationally symmetric in x-z."
    )


def test_positivity():
    """
    After running the Sod tube, density and pressure must be strictly positive
    everywhere in the domain.  Negative values would indicate a catastrophic
    numerical failure.
    """
    from tests.runner import get_pressure

    lo, hi, res, x_cell = make_slab_grid(N_SYM, SOD_T, interface_axis=0)
    h = 1.0 / N_SYM

    left = x_cell < SOD_X0
    rho_1d = np.where(left, SOD_RHO_L, SOD_RHO_R)
    p_1d = np.where(left, SOD_P_L, SOD_P_R)
    rho = rho_1d[:, None, None] * np.ones(res)
    p = p_1d[:, None, None] * np.ones(res)
    zeros = np.zeros(res)

    init_data = make_conserved(rho, zeros, zeros, zeros, p)
    _, _, _, data = run_test(
        name="SodPositivity",
        lo=lo, hi=hi, res=res,
        init_data=init_data,
        final_time=SOD_T,
        bc_lo=[BC_TRANSMISSIVE, BC_PERIODIC, BC_PERIODIC],
        bc_hi=[BC_TRANSMISSIVE, BC_PERIODIC, BC_PERIODIC],
    )

    rho_final = data[..., 0]
    p_final = get_pressure(data)

    assert np.all(rho_final > 0), f"Non-positive density found: min = {rho_final.min():.4e}"
    assert np.all(p_final > 0), f"Non-positive pressure found: min = {p_final.min():.4e}"
