from matplotlib import pyplot as plt
import numpy as np
from pathlib import Path

from file_handler import *
from mesh import *

test_out_dir = "test-output/"


def get_sdf_fname(fname):
    path = Path(fname)
    return str(path.parent / f"{path.stem}SDF{path.suffix}")


def plot_sod_test():
    name = "SodTest"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    rho = data[..., 0]

    plt.figure(figsize=(8, 6))
    plt.xlim(lo[0], hi[0])
    plt.plot(x, rho, ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".pdf")
    plt.close()


def plot_cylindrical_explosion():
    name = "CylindricalExplosion"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = data[..., 0]
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.xlim(lo[0], hi[0])
    plt.plot(x, rho[..., 0], ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".pdf")
    plt.close()


def plot_spherical_explosion():
    name = "SphericalExplosion"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    z = np.linspace(lo[2], hi[2], data.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = data[..., 0]

    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X[..., 0], Y[..., 0], rho[..., 0])
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(x, rho[..., 0, 0], ".")
    plt.xlim(lo[0], hi[0])
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".pdf")


def plot_kelvin_helmholtz():
    name = "KelvinHelmholtz"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = data[..., 0]
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho, cmap="Blues")
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()


def plot_shock_reflection():
    name = "ShockReflection"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    lo_sdf, hi_sdf, sdf = load_sdf(get_sdf_fname(remove_step_counter(get_last_header_fname(name, test_out_dir))))
    assert(np.allclose(lo_sdf, lo) and np.allclose(hi_sdf, hi))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = np.where(sdf[1:-1, 1:-1] < 0, np.nan, data[..., 0])
    
    plt.figure(figsize=(10, 5))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.contour(X, Y, sdf[1:-1, 1:-1], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()


def plot_hypersonic_sphere(useSTL):
    name = "HypersonicSphere"
    if useSTL:
        name += "FromSTL"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    lo_sdf, hi_sdf, sdf = load_sdf(get_sdf_fname(remove_step_counter(get_last_header_fname(name, test_out_dir))))
    assert(np.allclose(lo_sdf, lo) and np.allclose(hi_sdf, hi))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    z = np.linspace(lo[2], hi[2], data.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = np.where(sdf[1:-1, 1:-1, 1:-1] < 0, np.nan, data[..., 0])
    
    sliceIdx = data.shape[2] // 2
    plt.figure(figsize=(5, 6))
    plt.pcolormesh(X[..., sliceIdx], Y[..., sliceIdx], rho[..., sliceIdx])
    plt.colorbar(label="Density")
    plt.contour(X[..., sliceIdx], Y[..., sliceIdx], sdf[1:-1, 1:-1, sliceIdx], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()


def plot_space_shuttle():
    name = "SpaceShuttle"
    step, time, lo, hi, data = load(get_last_header_fname(name, test_out_dir))
    lo_sdf, hi_sdf, sdf = load_sdf(get_sdf_fname(remove_step_counter(get_last_header_fname(name, test_out_dir))))
    assert(np.allclose(lo_sdf, lo) and np.allclose(hi_sdf, hi))
    x = np.linspace(lo[0], hi[0], data.shape[0])
    y = np.linspace(lo[1], hi[1], data.shape[1])
    z = np.linspace(lo[2], hi[2], data.shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = np.where(sdf[1:-1, 1:-1, 1:-1] < 0, np.nan, data[..., 0])
    
    sliceIdx = data.shape[2] // 2
    plt.figure(figsize=(12, 10))
    plt.pcolormesh(X[..., sliceIdx], Y[..., sliceIdx], rho[..., sliceIdx])
    plt.colorbar(label="Density")
    plt.contour(X[..., sliceIdx], Y[..., sliceIdx], sdf[1:-1, 1:-1, sliceIdx], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(name + " at time " + str(time) + " after " + str(step) + " steps")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(test_out_dir + name + ".png", dpi=300)
    plt.close()


if __name__ == "__main__":
    try:
        plot_sod_test()
    except FileNotFoundError:
        pass

    try:
        plot_cylindrical_explosion()
    except FileNotFoundError:
        pass

    try:
        plot_spherical_explosion()
    except FileNotFoundError:
        pass

    try:
        plot_kelvin_helmholtz()
    except FileNotFoundError:
        pass

    try:
        plot_shock_reflection()
    except FileNotFoundError:
        pass

    try:
        plot_hypersonic_sphere(False)
    except FileNotFoundError:
        pass

    try:
        plot_hypersonic_sphere(True)
    except FileNotFoundError:
        pass

    try:
        plot_space_shuttle()
    except FileNotFoundError:
        pass
