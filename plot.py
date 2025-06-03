from matplotlib import pyplot as plt
import numpy as np


def read_output(fname):
    info = {}
    with open(fname, "r") as f:
        info['step'] = int(f.readline())
        info['time'] = float(f.readline())
        info['lo'] = [float(i) for i in f.readline().split()]
        info['hi'] = [float(i) for i in f.readline().split()]
        info['res'] = [int(i) for i in f.readline().split()]
        info['NVARS'] = int(f.readline())
        info['GRIDDIM'] = len(info['lo'])
        assert info['GRIDDIM'] == len(info['hi']) == len(info['res'])
        info['SPACEDIM'] = info['NVARS'] - 2
    
    return np.loadtxt(fname, skiprows=6).reshape((*info['res'], info['NVARS'])), info


def plot_sod_test():
    name = "SodTest"
    data, info = read_output(name + ".txt")
    x = np.linspace(info['lo'][0], info['hi'][0], info['res'][0])
    rho = data[:, 0]

    plt.figure(figsize=(8, 6))
    plt.xlim(info['lo'][0], info['hi'][0])
    plt.plot(x, rho, ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + ".pdf")
    plt.close()


def plot_cylindrical_explosion():
    name = "CylindricalExplosion"
    data, info = read_output(name + ".txt")
    x = np.linspace(info['lo'][0], info['hi'][0], info['res'][0])
    y = np.linspace(info['lo'][1], info['hi'][1], info['res'][1])
    X, Y = np.meshgrid(x, y, indexing='ij')
    rho = data[:, :, 0].transpose()
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.xlim(info['lo'][0], info['hi'][0])
    plt.plot(x, rho[:, 0], ".")
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + "Slice.pdf")
    plt.close()


def plot_spherical_explosion():
    name = "SphericalExplosion"
    data, info = read_output(name + ".txt")
    x = np.linspace(info['lo'][0], info['hi'][0], info['res'][0])
    y = np.linspace(info['lo'][1], info['hi'][1], info['res'][1])
    z = np.linspace(info['lo'][2], info['hi'][2], info['res'][2])
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    rho = data[:, :, :, 0].transpose((2, 1, 0))

    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X[:, :, 0], Y[:, :, 0], rho[:, :, 0])
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(name + ".png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(x, rho[:, 0, 0], ".")
    plt.xlim(info['lo'][0], info['hi'][0])
    plt.xlabel("x")
    plt.ylabel("Density")
    plt.tight_layout()
    plt.savefig(name + "Slice.pdf")


def plot_kelvin_helmholtz():
    name = "KelvinHelmholtz"
    data, info = read_output(name + ".txt")
    x = np.linspace(info['lo'][0], info['hi'][0], info['res'][0])
    y = np.linspace(info['lo'][1], info['hi'][1], info['res'][1])
    X, Y = np.meshgrid(x, y, indexing='ij')
    rho = data[:, :, 0].transpose()
    
    plt.figure(figsize=(8, 6))
    plt.pcolormesh(X, Y, rho)
    plt.colorbar(label="Density")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.gca().set_aspect('equal', adjustable='box')
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
