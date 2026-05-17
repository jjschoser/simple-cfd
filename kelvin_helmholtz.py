from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np
import subprocess

from mesh import *
from settings import *

# In order for this script to work, we require REAL double, 
# GRIDDIM 2, and SPACEDIM 2 in Macros.H

if __name__ == "__main__":
    test_out_dir = "test-output/kelvin-helmholtz-anim/"
    name = "KelvinHelmholtz"
    settings_fname = name + "Settings.txt"
    init_header_fname = name + ".txt"

    final_time = 3.0
    movie_duration = 10.0
    movie_fps = 25
    out_interval = final_time / (movie_duration * movie_fps)

    gamma = 1.4
    band_height = 0.5

    bc_lo = [2, 2]
    bc_hi = [2, 2]

    lo = np.array([-0.5, -0.5])
    hi = np.array([0.5, 0.5])
    res = np.array([512, 512])

    NVARS = 4
    init_data = np.zeros((*res, NVARS))
    dx = (hi - lo) / res
    x, y = [np.linspace(lo[d] + 0.5 * dx[d], hi[d] - 0.5 * dx[d], res[d]) for d in range(2)]
    X, Y = np.meshgrid(x, y, indexing="ij")

    init_data[:, :, 0] = np.where(np.abs(Y) < 0.5 * band_height, 2.0, 1.0)
    vel_x = np.where(np.abs(Y) < 0.5 * band_height, -0.5, 0.5)
    vel_y = 0.01 * np.sin(2.0 * np.pi * X)
    p = np.full_like(X, 2.5)
    init_data[:, :, 1] = init_data[:, :, 0] * vel_x
    init_data[:, :, 2] = init_data[:, :, 0] * vel_y
    init_data[:, :, 3] = 0.5 * init_data[:, :, 0] * (vel_x ** 2 + vel_y ** 2) + p / (gamma - 1)
    save(test_out_dir + init_header_fname, 0, 0.0, lo, hi, init_data)

    write_settings(test_out_dir + settings_fname, init_header_fname, 
                   final_time, bc_lo=bc_lo, bc_hi=bc_hi, gamma=gamma, 
                   out_interval=out_interval)
    subprocess.run(["make", "clean"])
    subprocess.run(["make"])
    subprocess.run(["./ilmatar-cfd", test_out_dir + settings_fname])

    header_fname_list = get_header_fname_list(name, test_out_dir)

    def load_frame_info(frame):
        step, time, _, __, data = load(header_fname_list[frame])
        return step, time, data[..., 0]
    
    rho_min = 0.5
    rho_max = 2.2
    step, time, rho = load_frame_info(0)
    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = plt.pcolormesh(X, Y, rho, cmap="Blues", vmin=rho_min, vmax=rho_max)
    plt.colorbar(mesh, label="Density")
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
