"""Smoke tests for the surface enumeration workflow."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
ENUM_INPUT = (
    PROJECT_ROOT
    / "examples"
    / "enumeration-examples"
    / "input"
    / "input-GaAs.json"
)
SYMMETRIC_ENUM_INPUT = (
    PROJECT_ROOT
    / "examples"
    / "enumeration-examples"
    / "input"
    / "input-Li2O.json"
)
MAKESTR_CANDIDATES = ("makestr.x", "makeStr.x", "makeStr.py")


def _enumlib_env(tmp_path: Path):
    """Return executable paths and an environment that can find enumlib."""
    env = os.environ.copy()
    executable_dir = Path(sys.executable).resolve().parent
    env["PATH"] = f"{executable_dir}{os.pathsep}{env.get('PATH', '')}"
    env["MPLCONFIGDIR"] = str(tmp_path / "mplconfig")

    enum_cmd = shutil.which("enum.x", path=env["PATH"])
    makestr_cmd = None
    for candidate in MAKESTR_CANDIDATES:
        makestr_cmd = shutil.which(candidate, path=env["PATH"])
        if makestr_cmd is not None:
            break
    return enum_cmd, makestr_cmd, env


def _skip_without_enumlib(tmp_path: Path):
    """Skip tests that require enumlib when the executables are unavailable."""
    enum_cmd, makestr_cmd, env = _enumlib_env(tmp_path)
    if enum_cmd is None or makestr_cmd is None:
        pytest.skip(
            "enumlib executables are required: expected enum.x and one of "
            f"{', '.join(MAKESTR_CANDIDATES)}"
        )
    return enum_cmd, makestr_cmd, env


def _write_absolute_input(source_path: Path, tmp_path: Path):
    """Write a temp enumeration input with absolute structure paths."""
    data = json.loads(source_path.read_text())
    data["target_slab_path"] = str(
        (PROJECT_ROOT / data["target_slab_path"]).resolve()
    )
    input_path = tmp_path / source_path.name
    input_path.write_text(json.dumps(data))
    return input_path


def _run_surface_enumeration(input_path: Path, tmp_path: Path, env: dict):
    """Run the surface-enumeration module against an input file."""
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "surface_pd.cli.surface_enumeration",
            str(input_path),
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )


def test_enumlib_executables_are_available(tmp_path):
    """The enumlib executable pair should be discoverable when installed."""
    enum_cmd, makestr_cmd, _ = _skip_without_enumlib(tmp_path)

    assert Path(enum_cmd).is_file()
    assert Path(makestr_cmd).is_file()


def test_surface_enumeration_cli_smoke(tmp_path):
    """Run a small deterministic enumeration through the CLI entry path."""
    _, _, env = _skip_without_enumlib(tmp_path)
    input_path = _write_absolute_input(ENUM_INPUT, tmp_path)

    result = _run_surface_enumeration(input_path, tmp_path, env)

    assert result.returncode == 0, result.stderr
    assert "target_cell_size = 2" in result.stdout
    assert (
        "The enumeration found 20(20+0) distinct structures for "
        "['Ga', 'As'] with [0.75, 0.75] composition."
    ) in result.stdout
    assert "20 distinct structures are found totally." in result.stdout


def test_symmetric_surface_enumeration_cli_smoke(tmp_path):
    """Run a deterministic symmetric enumeration through the CLI path."""
    _, _, env = _skip_without_enumlib(tmp_path)
    input_path = _write_absolute_input(SYMMETRIC_ENUM_INPUT, tmp_path)

    result = _run_surface_enumeration(input_path, tmp_path, env)

    assert result.returncode == 0, result.stderr
    assert "target_cell_size = 2" in result.stdout
    assert (
        "The enumeration found 6(6+0) distinct structures for "
        "['Li', 'O'] with [1.0, 0.75] composition."
    ) in result.stdout
    assert "6 distinct structures are found totally." in result.stdout
