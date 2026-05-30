"""Build and execution utilities for ilmatar-cfd tests."""

import re
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
MACROS_H = PROJECT_ROOT / "Macros.H"
BINARY = PROJECT_ROOT / "ilmatar-cfd"
TEST_OUT_DIR = PROJECT_ROOT / "test-output" / "pytest"

sys.path.insert(0, str(PROJECT_ROOT))

from mesh import save, load  # noqa: E402
from file_handler import get_last_header_fname  # noqa: E402
from settings import write_settings  # noqa: E402

BC_TRANSMISSIVE = 0
BC_REFLECTIVE = 1
BC_PERIODIC = 2

DEFAULT_GAMMA = 1.4


def read_macros() -> dict:
    content = MACROS_H.read_text()
    griddim = int(re.search(r"^#define GRIDDIM (\d+)", content, re.MULTILINE).group(1))
    spacedim = int(re.search(r"^#define SPACEDIM (\d+)", content, re.MULTILINE).group(1))
    use_rigid = bool(re.search(r"^\s*#define USE_RIGID\b", content, re.MULTILINE))
    return {"griddim": griddim, "spacedim": spacedim, "use_rigid": use_rigid}


def make_conserved(rho, ux, uy, uz, p, gamma=DEFAULT_GAMMA):
    """Primitive → conserved for GRIDDIM=3, SPACEDIM=3. Last axis is NVARS=5."""
    rhoE = p / (gamma - 1) + 0.5 * rho * (ux**2 + uy**2 + uz**2)
    return np.stack([rho, rho * ux, rho * uy, rho * uz, rhoE], axis=-1)


def get_pressure(data, gamma=DEFAULT_GAMMA):
    """Recover pressure from conserved variable array (last axis = NVARS=5)."""
    rho = data[..., 0]
    rhoE = data[..., 4]
    kin = 0.5 * (data[..., 1] ** 2 + data[..., 2] ** 2 + data[..., 3] ** 2) / rho
    return (gamma - 1) * (rhoE - kin)


def run_test(
    name, lo, hi, res, init_data, final_time, bc_lo, bc_hi,
    gamma=DEFAULT_GAMMA, cfl=None,
):
    """
    Write initial data and settings, run the solver, return (time, lo, hi, data).
    All output files land in TEST_OUT_DIR.
    name must be unique across concurrent tests to avoid file collisions.
    """
    TEST_OUT_DIR.mkdir(parents=True, exist_ok=True)

    init_fname = name + "Init.txt"
    settings_fname = name + "Settings.txt"
    out_base_fname = name + ".txt"

    save(str(TEST_OUT_DIR / init_fname), 0, 0.0, lo, hi, init_data)

    write_settings(
        str(TEST_OUT_DIR / settings_fname),
        init_header_fname=init_fname,
        final_time=final_time,
        bc_lo=bc_lo,
        bc_hi=bc_hi,
        gamma=gamma,
        cfl=cfl,
        out_header_base_fname=out_base_fname,
    )

    result = subprocess.run(
        [str(BINARY), str(TEST_OUT_DIR / settings_fname)],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Solver failed for '{name}':\n{result.stderr}")

    _, time, lo_out, hi_out, data = load(
        get_last_header_fname(name, str(TEST_OUT_DIR))
    )
    return time, lo_out, hi_out, data


def make_slab_grid(N, final_time=None, interface_axis=0):
    """
    Build a thin 3D slab suited for 1D-in-3D tests.
    Returns (lo, hi, res, x_cell) where x_cell are cell centres along interface_axis.

    The slab has N cells along the chosen axis and 4 cells along each transverse
    axis.  All cells are square: dx = dy = dz = 1/N.
    """
    h = 1.0 / N
    transverse = 4 * h

    lo = [0.0, 0.0, 0.0]
    hi = [transverse, transverse, transverse]
    res = [4, 4, 4]

    hi[interface_axis] = 1.0
    res[interface_axis] = N

    x_cell = np.array([lo[interface_axis] + (i + 0.5) * h for i in range(N)])
    return lo, hi, res, x_cell
