from matplotlib import pyplot as plt
import numpy as np
import subprocess

from mesh import *
from settings import *

# In order for this script to work, we require REAL double, 
# GRIDDIM 2, SPACEDIM 2, and USE_RIGID in Macros.H

if __name__ == "__main__":
    test_out_dir = "test-output/"
    name = "ShockReflectionFromScript"
    settings_fname = name + "Settings.txt"
    init_header_fname = name + ".txt"
    sdf_header_fname = name + "SDF.txt"

    gamma = 1.4
    x_shock = 4e-3
    x_wedge = 4.96e-3
    alpha_wedge = np.deg2rad(25.0)
    final_time = 35e-6
    rho_inf = 1.225
    vel_inf = 0
    p_inf = 101325
    M_shock = 1.7

    bc_lo = [0, 1]
    bc_hi = [0, 0]

    lo = np.array([-4e-3, 0.0])
    hi = np.array([29e-3, 16.5e-3])
    res = np.array([1024, 512])
    
    c_inf = np.sqrt(gamma * p_inf / rho_inf)
    s_shock = M_shock * c_inf
    M_inf = vel_inf / c_inf
    rho_star = rho_inf * ((gamma + 1) * (M_inf - M_shock) ** 2) / ((gamma - 1) * (M_inf - M_shock) ** 2 + 2)
    p_star = p_inf * ((2 * gamma * (M_inf - M_shock) ** 2 - (gamma - 1)) / (gamma + 1))
    vel_star = (1 - rho_inf / rho_star) * s_shock + vel_inf * (1 - rho_inf / rho_star)

    NVARS = 4
    init_data = np.zeros((*res, NVARS))
    sdf = np.zeros(res)
    dx = (hi - lo) / res
    x, y = [np.linspace(lo[d] + 0.5 * dx[d], hi[d] - 0.5 * dx[d], res[d]) for d in range(2)]
    X, Y = np.meshgrid(x, y, indexing="ij")

    init_data[:, :, 0] = np.where(X < x_shock, rho_star, rho_inf)
    vel_x = np.where(X < x_shock, vel_star, vel_inf)
    vel_y = np.zeros_like(X)
    p = np.where(X < x_shock, p_star, p_inf)
    init_data[:, :, 1] = init_data[:, :, 0] * vel_x
    init_data[:, :, 2] = init_data[:, :, 0] * vel_y
    init_data[:, :, 3] = 0.5 * init_data[:, :, 0] * (vel_x ** 2 + vel_y ** 2) + p / (gamma - 1)
    save(test_out_dir + init_header_fname, 0, 0.0, lo, hi, init_data)

    x_sdf, y_sdf = [np.linspace(lo[d] - 0.5 * dx[d], hi[d] + 0.5 * dx[d], res[d] + 2) for d in range(2)]
    X_sdf, Y_sdf = np.meshgrid(x_sdf, y_sdf, indexing="ij")
    sdf = np.cos(alpha_wedge) * Y_sdf - np.sin(alpha_wedge) * (X_sdf - x_wedge)
    save_sdf(test_out_dir + sdf_header_fname, lo, hi, sdf)

    write_settings(test_out_dir + settings_fname, init_header_fname, 
                   final_time, bc_lo=bc_lo, bc_hi=bc_hi, gamma=gamma, 
                   sdf_header_fname=sdf_header_fname)
    subprocess.run(["make", "clean"])
    subprocess.run(["make"])
    subprocess.run(["./ilmatar-cfd", test_out_dir + settings_fname])

    step, time, _, __, final_data = load(get_last_header_fname(name, test_out_dir))
    rho = final_data[:, :, 0]
    schlieren = np.where(sdf[1:-1, 1:-1] < 0, 0.0, np.exp(-np.linalg.norm(np.asarray(np.gradient(rho, *dx)), axis=0)/np.sqrt(rho)/2000))
    
    plt.figure(figsize=(10, 5))
    plt.pcolormesh(X, Y, schlieren, cmap="gray")
    plt.contour(X, Y, sdf[1:-1, 1:-1], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()
