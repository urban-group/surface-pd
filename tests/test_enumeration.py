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
LI_SLAB = (
    PROJECT_ROOT
    / "examples"
    / "enumeration-examples"
    / "structure"
    / "electrode"
    / "POSCAR_Li_100.vasp"
)
LI2O_SLAB = (
    PROJECT_ROOT
    / "examples"
    / "enumeration-examples"
    / "structure"
    / "electrode"
    / "POSCAR_Li2O_110.vasp"
)


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
        (PROJECT_ROOT / "examples" / data["target_slab_path"]).resolve()
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


def test_public_enumerator_returns_only_in_plane_surface_cells(tmp_path):
    """Reject enumlib derivatives that multiply or tilt the slab normal."""
    _, _, env = _skip_without_enumlib(tmp_path)
    script = f"""
from pymatgen.core import Structure
from surface_pd.core import EnumerationSlab, SurfaceEnumerator

source = Structure.from_file({str(LI_SLAB)!r})
slab = EnumerationSlab.from_structure(
    source,
    direction=2,
    enumerated_species=["Li"],
    num_enumerated_layers={{"Li": 1}},
    symmetric=False,
)
results = SurfaceEnumerator(
    {{"Li": {{"Li": 0.5}}}}, min_cell_size=2, max_cell_size=2
).apply_enumeration(slab, max_structures=20)
print(len(results))
for result in results:
    print(len(result), *result.lattice.abc)
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )

    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0] == "2"
    structures = [line.split() for line in lines[1:]]
    assert [int(values[0]) for values in structures] == [17, 17]
    assert all(
        float(values[3]) == pytest.approx(30.791121)
        for values in structures
    )


def test_public_enumerator_finalizes_symmetric_surface_cells(tmp_path):
    """Return complete inversion-symmetric Li2O surface derivatives."""
    _, _, env = _skip_without_enumlib(tmp_path)
    script = f"""
from pymatgen.core import Structure
from surface_pd.core import EnumerationSlab, SurfaceEnumerator

source = Structure.from_file({str(LI2O_SLAB)!r})
slab = EnumerationSlab.from_structure(
    source,
    direction=2,
    enumerated_species=["Li", "O"],
    num_enumerated_layers={{"Li": 1, "O": 2}},
    symmetric=True,
)
results = SurfaceEnumerator(
    {{"Li": {{"Li": 1.0}}, "O": {{"O": 0.75}}}},
    min_cell_size=1,
    max_cell_size=2,
).apply_enumeration(slab, max_structures=20)
print(len(results))
for result in results:
    complete = all(
        site.properties.get("selective_dynamics") is not None
        for site in result
    )
    print(
        len(result), result.composition["Li"], result.composition["O"],
        result.lattice.c, result.has_inversion_symmetry(), complete,
    )
"""
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )

    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0] == "4"
    for line in lines[1:]:
        sites, lithium, oxygen, c_length, symmetric, complete = line.split()
        assert (int(sites), float(lithium), float(oxygen)) == (76, 52.0, 24.0)
        assert float(c_length) == pytest.approx(38.7623252869)
        assert symmetric == "True"
        assert complete == "True"
