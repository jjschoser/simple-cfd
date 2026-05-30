"""
Conservation tests.

For periodic boundary conditions, the Euler finite volume scheme conserves
mass, momentum, and energy exactly (to floating-point precision) via telescoping
cancellation of fluxes.  These tests verify that property.
"""

import numpy as np
import pytest

from tests.runner import (
    BC_PERIODIC,
    DEFAULT_GAMMA,
    make_conserved,
    run_test,
)


def _initial_data():
    """
    16³ cell periodic domain with a smooth density wave and non-zero background
    velocity.  The non-zero velocity ensures mass flux crosses cell faces,
    exercising the flux computation in all three directions.
    """
    gamma = DEFAULT_GAMMA
    N = 16
    lo = [0.0, 0.0, 0.0]
    hi = [1.0, 1.0, 1.0]
    res = [N, N, N]

    h = 1.0 / N
    x = lo[0] + (np.arange(N) + 0.5) * h
    y = lo[1] + (np.arange(N) + 0.5) * h
    z = lo[2] + (np.arange(N) + 0.5) * h
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    rho = 1.0 + 0.1 * np.sin(2 * np.pi * X)
    ux = 0.3 * np.ones_like(X)
    uy = 0.1 * np.ones_like(X)
    uz = np.zeros_like(X)
    p = np.ones_like(X)

    data = make_conserved(rho, ux, uy, uz, p, gamma)
    return lo, hi, res, data


def test_conservation_periodic():
    """
    After running in a periodic box, each conserved integral
    (mass, x-momentum, y-momentum, z-momentum, total energy)
    must match its initial value to within 1e-10 relative tolerance.
    """
    gamma = DEFAULT_GAMMA
    lo, hi, res, init_data = _initial_data()
    bc_all = [BC_PERIODIC, BC_PERIODIC, BC_PERIODIC]

    # Record initial sums (proportional to cell-volume integrals, same dx³ everywhere)
    sum_init = init_data.sum(axis=(0, 1, 2))  # shape (5,)

    _, _, _, final_data = run_test(
        name="Conservation",
        lo=lo,
        hi=hi,
        res=res,
        init_data=init_data,
        final_time=0.5,
        bc_lo=bc_all,
        bc_hi=bc_all,
        gamma=gamma,
    )

    sum_final = final_data.sum(axis=(0, 1, 2))

    labels = ["mass", "x-momentum", "y-momentum", "z-momentum", "total energy"]
    print("\nConservation errors (periodic box, t=0.5):")
    print(f"  {'Quantity':<15}  {'Initial sum':>14}  {'Final sum':>14}  {'Rel. error':>12}")
    print(f"  {'-'*15}  {'-'*14}  {'-'*14}  {'-'*12}")
    for i, label in enumerate(labels):
        rel_err = abs(sum_final[i] - sum_init[i]) / (abs(sum_init[i]) + 1e-30)
        print(f"  {label:<15}  {sum_init[i]:>14.6e}  {sum_final[i]:>14.6e}  {rel_err:>12.2e}")
        assert rel_err < 1e-10, (
            f"{label} not conserved: initial sum = {sum_init[i]:.6e}, "
            f"final sum = {sum_final[i]:.6e}, relative error = {rel_err:.2e}"
        )
