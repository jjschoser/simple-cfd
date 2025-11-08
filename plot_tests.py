from matplotlib import pyplot as plt
import numpy as np


def read_output(fname, double=True):
    info = {}
    with open(fname, "r") as f:
        info["step"] = int(f.readline())
        info["time"] = float(f.readline())
        info["lo"] = [float(i) for i in f.readline().split()]
        info["hi"] = [float(i) for i in f.readline().split()]
        info["res"] = [int(i) for i in f.readline().split()]
        info["NVARS"] = int(f.readline())
        info["data_filename"] = f.readline().rstrip()
        info["GRIDDIM"] = len(info["lo"])
        assert info["GRIDDIM"] == len(info["hi"]) == len(info["res"])
        info["SPACEDIM"] = info["NVARS"] - 2

    REAL = np.float64 if double else np.float32
    count = np.prod(info["res"]) * info["NVARS"]
    with open(info["data_filename"], "rb") as f:
        data = np.fromfile(f, dtype=REAL, count=count).reshape((*info["res"], info["NVARS"]))

    try:
        sdf_count = np.prod(np.asarray(info["res"]) + 2)
        with open(info["data_filename"].replace(".dat", "SDF.dat"), "rb") as f:
            sdf_data = np.fromfile(f, dtype=REAL, count=sdf_count).reshape(np.asarray(info["res"]) + 2)
        return data, sdf_data, info
    except FileNotFoundError:
        return data, info


def plot_sod_test():
    name = "SodTest"
    data, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    rho = data[..., 0]

    plt.figure(figsize=(8, 6))
    plt.xlim(info["lo"][0], info["hi"][0])
    plt.plot(x, rho, ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + ".pdf")
    plt.close()


def plot_cylindrical_explosion():
    name = "CylindricalExplosion"
    data, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = data[..., 0]
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.xlim(info["lo"][0], info["hi"][0])
    plt.plot(x, rho[..., 0], ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + "Slice.pdf")
    plt.close()


def plot_spherical_explosion():
    name = "SphericalExplosion"
    data, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    z = np.linspace(info["lo"][2], info["hi"][2], info["res"][2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = data[..., 0]

    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X[..., 0], Y[..., 0], rho[..., 0])
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(x, rho[..., 0, 0], ".")
    plt.xlim(info["lo"][0], info["hi"][0])
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + "Slice.pdf")


def plot_kelvin_helmholtz():
    name = "KelvinHelmholtz"
    data, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = data[..., 0]
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho, cmap="Blues")
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()


def plot_shock_reflection():
    name = "ShockReflection"
    data, sdf, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    rho = np.where(sdf[1:-1, 1:-1] < 0, np.nan, data[..., 0])
    
    plt.figure(figsize=(10, 5))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.contour(X, Y, sdf[1:-1, 1:-1], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()


def plot_hypersonic_sphere(useSTL):
    name = "HypersonicSphere"
    if useSTL:
        name += "FromSTL"
    data, sdf, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    z = np.linspace(info["lo"][2], info["hi"][2], info["res"][2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = np.where(sdf[1:-1, 1:-1, 1:-1] < 0, np.nan, data[..., 0])
    
    sliceIdx = info["res"][2] // 2
    plt.figure(figsize=(5, 6))
    plt.pcolormesh(X[..., sliceIdx], Y[..., sliceIdx], rho[..., sliceIdx])
    plt.colorbar(label="Density")
    plt.contour(X[..., sliceIdx], Y[..., sliceIdx], sdf[1:-1, 1:-1, sliceIdx], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()


def plot_wing():
    name = "Wing"
    data, sdf, info = read_output(name + ".txt")
    x = np.linspace(info["lo"][0], info["hi"][0], info["res"][0])
    y = np.linspace(info["lo"][1], info["hi"][1], info["res"][1])
    z = np.linspace(info["lo"][2], info["hi"][2], info["res"][2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    rho = data[..., 0]
    mom = data[..., 1:4]
    velz = np.where(sdf[1:-1, 1:-1, 1:-1] < 0, np.nan, mom[...,2]/rho)
    velzmax = np.max(np.abs(velz))
    
    sliceIdx = info["res"][2] // 2
    plt.figure(figsize=(12, 6))
    plt.pcolormesh(X[..., sliceIdx], Y[..., sliceIdx], velz[..., sliceIdx], cmap="coolwarm", vmin=-velzmax, vmax=velzmax)
    plt.colorbar(label="z-velocity")
    plt.contour(X[..., sliceIdx], Y[..., sliceIdx], sdf[1:-1, 1:-1, sliceIdx], levels=[0], colors="k")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
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
        plot_wing()
    except FileNotFoundError:
        pass
