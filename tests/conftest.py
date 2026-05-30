"""Shared pytest fixtures for ilmatar-cfd tests."""

import re
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
MACROS_H = PROJECT_ROOT / "Macros.H"
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture(scope="session", autouse=True)
def ensure_griddim3_build():
    """
    Guarantee the solver binary is compiled with GRIDDIM=3, SPACEDIM=3.
    If the current Macros.H differs, patch it, rebuild, and restore on teardown.
    All tests in this suite require the 3-D binary.
    """
    from tests.runner import read_macros

    original_content = MACROS_H.read_text()
    macros = read_macros()

    needs_rebuild = macros["griddim"] != 3 or macros["spacedim"] != 3

    if needs_rebuild:
        content = original_content
        content = re.sub(r"(#define GRIDDIM )\d+", r"\g<1>3", content)
        content = re.sub(r"(#define SPACEDIM )\d+", r"\g<1>3", content)
        MACROS_H.write_text(content)
        subprocess.run(["make", "clean"], cwd=PROJECT_ROOT, check=True, capture_output=True)
        subprocess.run(["make"], cwd=PROJECT_ROOT, check=True, capture_output=True)

    yield

    if needs_rebuild:
        MACROS_H.write_text(original_content)
        # Rebuild so the project is left in its original state.
        subprocess.run(["make", "clean"], cwd=PROJECT_ROOT, check=True, capture_output=True)
        subprocess.run(["make"], cwd=PROJECT_ROOT, check=True, capture_output=True)
