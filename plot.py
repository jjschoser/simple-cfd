from matplotlib import pyplot as plt
import numpy as np

if __name__ == "__main__":
    fname = "out.txt"

    with open(fname, "r") as f:
        step = int(f.readline())
        time = float(f.readline())
        lo = [float(i) for i in f.readline().split()]
        hi = [float(i) for i in f.readline().split()]
        res = [int(i) for i in f.readline().split()]
        NVARS = int(f.readline())
        GRIDDIM = len(lo)
        assert GRIDDIM == len(hi) == len(res)
        SPACEDIM = NVARS - 2
    
    data = np.loadtxt(fname, skiprows=6).reshape((*res, NVARS))

    if GRIDDIM == 1:
        x = np.linspace(lo[0], hi[0], res[0])
        rho = data[:, 0]

        plt.figure(figsize=(8, 6))
        plt.xlim(lo[0], hi[0])
        plt.plot(x, rho, ".")
        plt.xlabel("x")
        plt.ylabel("Density")
        plt.tight_layout()
        plt.savefig("out.pdf")

    elif GRIDDIM == 2:
        x = np.linspace(lo[0], hi[0], res[0])
        y = np.linspace(lo[1], hi[1], res[1])
        X, Y = np.meshgrid(x, y, indexing='ij')
        rho = data[:, :, 0]

        plt.figure(figsize=(8, 6))
        plt.pcolormesh(X, Y, rho)
        plt.colorbar(label="Density")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.savefig("out.png", dpi=300)
        plt.close()

        plt.figure(figsize=(8, 6))
        plt.xlim(lo[0], hi[0])
        plt.plot(x, rho[:, 0], ".")
        plt.xlabel("x")
        plt.ylabel("Density")
        plt.tight_layout()
        plt.savefig("out.pdf")
    
    else:
        x = np.linspace(lo[0], hi[0], res[0])
        y = np.linspace(lo[1], hi[1], res[1])
        z = np.linspace(lo[2], hi[2], res[2])
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        rho = data[:, :, :, 0]

        fig = plt.figure(figsize=(8, 6))
        plt.pcolormesh(X[:, :, 0], Y[:, :, 0], rho[:, :, 0])
        plt.colorbar(label="Density")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        plt.savefig("out.png", dpi=300)
        plt.close()

        fig = plt.figure(figsize=(8, 6))
        plt.plot(x, rho[:, 0, 0], ".")
        plt.xlim(lo[0], hi[0])
        plt.xlabel("x")
        plt.ylabel("Density")
        plt.tight_layout()
        plt.savefig("out.pdf")
