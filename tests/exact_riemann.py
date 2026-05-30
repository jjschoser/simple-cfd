"""
Exact solver for the 1D Riemann problem for the Euler equations.

Reference: Toro, E. F. (2009). Riemann Solvers and Numerical Methods for Fluid
           Dynamics (3rd ed.), Sections 4.2–4.4.
"""

import numpy as np


def _f(p, rho_k, p_k, a_k, gamma):
    """Wave function f_k(p): velocity change across the k-th wave (Toro 4.44)."""
    A = 2.0 / ((gamma + 1) * rho_k)
    B = (gamma - 1) / (gamma + 1) * p_k
    if p <= p_k:
        return 2 * a_k / (gamma - 1) * ((p / p_k) ** ((gamma - 1) / (2 * gamma)) - 1)
    else:
        return (p - p_k) * (A / (p + B)) ** 0.5


def _df(p, rho_k, p_k, a_k, gamma):
    """df_k/dp for Newton–Raphson."""
    A = 2.0 / ((gamma + 1) * rho_k)
    B = (gamma - 1) / (gamma + 1) * p_k
    if p <= p_k:
        return (p / p_k) ** (-(gamma + 1) / (2 * gamma)) / (rho_k * a_k)
    else:
        return (A / (p + B)) ** 0.5 * (1 - (p - p_k) / (2 * (p + B)))


def _solve_p_star(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma):
    """Find the contact pressure p* by Newton–Raphson (Toro 4.46)."""
    a_l = (gamma * p_l / rho_l) ** 0.5
    a_r = (gamma * p_r / rho_r) ** 0.5

    # Two-rarefaction initial guess
    nu = (gamma - 1) / (2 * gamma)
    p_star = (
        (a_l + a_r - (gamma - 1) / 2 * (u_r - u_l))
        / (a_l / p_l**nu + a_r / p_r**nu)
    ) ** (1.0 / nu)
    p_star = max(float(p_star), 1e-12)

    for _ in range(50):
        f_val = _f(p_star, rho_l, p_l, a_l, gamma) + _f(p_star, rho_r, p_r, a_r, gamma) + (u_r - u_l)
        df_val = _df(p_star, rho_l, p_l, a_l, gamma) + _df(p_star, rho_r, p_r, a_r, gamma)
        p_new = max(p_star - f_val / df_val, 1e-12)
        if abs(p_new - p_star) / (p_star + 1e-15) < 1e-10:
            return p_new
        p_star = p_new

    return p_star


def _sample_point(xi, rho_l, u_l, p_l, rho_r, u_r, p_r, p_star, u_star, gamma):
    """Sample the exact Riemann solution at a single xi = (x - x0) / t."""
    a_l = (gamma * p_l / rho_l) ** 0.5
    a_r = (gamma * p_r / rho_r) ** 0.5

    if xi <= u_star:
        # ── Left of contact ──────────────────────────────────────────────────
        if p_star <= p_l:
            # Left rarefaction
            a_l_star = a_l * (p_star / p_l) ** ((gamma - 1) / (2 * gamma))
            SHL = u_l - a_l
            STL = u_star - a_l_star
            if xi <= SHL:
                return rho_l, u_l, p_l
            elif xi <= STL:
                D = 2 / (gamma + 1) + (gamma - 1) / ((gamma + 1) * a_l) * (u_l - xi)
                rho = rho_l * D ** (2 / (gamma - 1))
                u = 2 / (gamma + 1) * (a_l + (gamma - 1) / 2 * u_l + xi)
                p = p_l * D ** (2 * gamma / (gamma - 1))
                return rho, u, p
            else:
                rho_star_l = rho_l * (p_star / p_l) ** (1 / gamma)
                return rho_star_l, u_star, p_star
        else:
            # Left shock
            SL = u_l - a_l * ((gamma + 1) / (2 * gamma) * p_star / p_l + (gamma - 1) / (2 * gamma)) ** 0.5
            if xi <= SL:
                return rho_l, u_l, p_l
            else:
                rho_star_l = rho_l * (
                    (p_star / p_l + (gamma - 1) / (gamma + 1))
                    / (1 + (gamma - 1) / (gamma + 1) * p_star / p_l)
                )
                return rho_star_l, u_star, p_star
    else:
        # ── Right of contact ─────────────────────────────────────────────────
        if p_star <= p_r:
            # Right rarefaction
            a_r_star = a_r * (p_star / p_r) ** ((gamma - 1) / (2 * gamma))
            SHR = u_r + a_r
            STR = u_star + a_r_star
            if xi >= SHR:
                return rho_r, u_r, p_r
            elif xi >= STR:
                D = 2 / (gamma + 1) - (gamma - 1) / ((gamma + 1) * a_r) * (u_r - xi)
                rho = rho_r * D ** (2 / (gamma - 1))
                u = 2 / (gamma + 1) * (-a_r + (gamma - 1) / 2 * u_r + xi)
                p = p_r * D ** (2 * gamma / (gamma - 1))
                return rho, u, p
            else:
                rho_star_r = rho_r * (p_star / p_r) ** (1 / gamma)
                return rho_star_r, u_star, p_star
        else:
            # Right shock
            SR = u_r + a_r * ((gamma + 1) / (2 * gamma) * p_star / p_r + (gamma - 1) / (2 * gamma)) ** 0.5
            if xi >= SR:
                return rho_r, u_r, p_r
            else:
                rho_star_r = rho_r * (
                    (p_star / p_r + (gamma - 1) / (gamma + 1))
                    / (1 + (gamma - 1) / (gamma + 1) * p_star / p_r)
                )
                return rho_star_r, u_star, p_star


def exact_riemann_1d(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma, x, t, x0=0.0):
    """
    Sample the exact Riemann solution at positions x and time t.

    Parameters
    ----------
    rho_l, u_l, p_l : float
        Left primitive state (density, velocity, pressure).
    rho_r, u_r, p_r : float
        Right primitive state.
    gamma : float
        Ratio of specific heats.
    x : array-like, shape (N,)
        Sample positions.
    t : float
        Sample time (must be > 0).
    x0 : float
        Position of the initial discontinuity.

    Returns
    -------
    rho, u, p : ndarray, shape (N,)
        Primitive variables at (x, t).
    """
    x = np.asarray(x, dtype=float)
    p_star = _solve_p_star(rho_l, u_l, p_l, rho_r, u_r, p_r, gamma)
    a_l = (gamma * p_l / rho_l) ** 0.5
    u_star = u_l - _f(p_star, rho_l, p_l, a_l, gamma)

    xi = (x - x0) / t
    rho = np.empty_like(x)
    u = np.empty_like(x)
    p = np.empty_like(x)

    for i, xi_i in enumerate(xi):
        rho[i], u[i], p[i] = _sample_point(
            xi_i, rho_l, u_l, p_l, rho_r, u_r, p_r, p_star, u_star, gamma
        )

    return rho, u, p
