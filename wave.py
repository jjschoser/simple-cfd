from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np
import subprocess

from mesh import *
from settings import *

# In order for this script to work, we require REAL double, 
# GRIDDIM 2, and SPACEDIM 2, USE_RIGID and USE_GRAVITY in Macros.H

if __name__ == "__main__":
    test_out_dir = "test-output/wave-anim/"
    name = "Wave"
    settings_fname = name + "Settings.txt"
    init_header_fname = name + ".txt"
    sdf_header_fname = name + "SDF.txt"

    final_time = 50.0
    movie_duration = 10.0
    movie_fps = 25
    out_interval = final_time / (movie_duration * movie_fps)

    gamma = 1.4

    bc_lo = [1, 1]
    bc_hi = [1, 1]

    lo = np.array([0.0, 0.0])
    hi = np.array([64.0, 16.0])
    res = np.array([1024, 256])

    NVARS = 4
    init_data = np.zeros((*res, NVARS))
    dx = (hi - lo) / res
    x, y = [np.linspace(lo[d] + 0.5 * dx[d], hi[d] - 0.5 * dx[d], res[d]) for d in range(2)]
    X, Y = np.meshgrid(x, y, indexing="ij")

    sea_level = 4.0
    rho_air = 0.01
    rho_water = 1.0
    g = 0.5
    p_sea_level = 1.0
    alpha_sea_floor = np.deg2rad(5.0)

    def is_ocean(y):
        return y < sea_level
    
    def is_droplet(x, y):
        return np.hypot(x, y - (sea_level + 4.0)) < 3.0

    init_data[:, :, 0] = np.where(np.logical_or(is_ocean(Y), is_droplet(X, Y)), rho_water, rho_air)
    vel_y = np.where(is_droplet(X, Y), -0.5, 0.0)
    p = np.where(is_ocean(Y), p_sea_level + rho_water * g * (sea_level - Y), p_sea_level - rho_air * g * (Y - sea_level))
    init_data[:, :, 1] = 0.0
    init_data[:, :, 2] = init_data[:, :, 0] * vel_y
    init_data[:, :, 3] = 0.5 * init_data[:, :, 0] * (vel_y ** 2) + p / (gamma - 1)
    save(test_out_dir + init_header_fname, 0, 0.0, lo, hi, init_data)

    x_sdf, y_sdf = [np.linspace(lo[d] - 0.5 * dx[d], hi[d] + 0.5 * dx[d], res[d] + 2) for d in range(2)]
    X_sdf, Y_sdf = np.meshgrid(x_sdf, y_sdf, indexing="ij")
    sdf = np.cos(alpha_sea_floor) * Y_sdf - np.sin(alpha_sea_floor) * (X_sdf - alpha_sea_floor)
    save_sdf(test_out_dir + sdf_header_fname, lo, hi, sdf)

    write_settings(test_out_dir + settings_fname, init_header_fname, 
                   final_time, bc_lo=bc_lo, bc_hi=bc_hi, gamma=gamma, 
                   sdf_header_fname=sdf_header_fname, out_interval=out_interval,
                   g=[0, -g])
    subprocess.run(["make", "clean"])
    subprocess.run(["make"])
    subprocess.run(["./ilmatar-cfd", test_out_dir + settings_fname])

    header_fname_list = get_header_fname_list(name, test_out_dir)

    def load_frame_info(frame):
        step, time, _, __, data = load(header_fname_list[frame])
        rho = np.where(sdf[1:-1, 1:-1] < 0, np.nan, data[:, :, 0])
        return step, time, rho
    
    step, time, rho = load_frame_info(0)
    fig, ax = plt.subplots(figsize=(15, 4))
    plt.contour(X, Y, sdf[1:-1, 1:-1], levels=[0], colors="k")
    mesh = plt.pcolormesh(X, Y, rho, cmap="Blues", vmin=0.9*rho_air, vmax=rho_water*1.5)
    plt.colorbar(mesh, label="Density", fraction=0.05)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{name} at time {time:.2f} after {step} steps")
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    
    def update(frame):
        step, time, rho = load_frame_info(frame)
        mesh.set_array(rho)
        ax.set_title(f"{name} at time {time:.2f} after {step} steps")
    
    print("Rendering animation...")
    anim = animation.FuncAnimation(fig, update, frames=len(header_fname_list), interval=1000/movie_fps)
    writer = animation.FFMpegWriter(fps=movie_fps)
    anim.save(test_out_dir + name + ".mp4", writer=writer, dpi=300)
    plt.close()
