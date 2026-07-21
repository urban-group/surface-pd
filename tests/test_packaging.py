"""Tests for built-wheel and installed-package behavior."""

import os
import subprocess
import sys
import venv
from pathlib import Path
from zipfile import ZipFile

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONSOLE_SCRIPTS = (
    "surface-enumeration",
    "surface-pd-plot",
    "generate-discharge-pd",
)
PACKAGE_FILES = (
    "surface_pd/_version.py",
    "surface_pd/__init__.py",
    "surface_pd/analysis/__init__.py",
    "surface_pd/analysis/slab_analysis.py",
    "surface_pd/cli/__init__.py",
    "surface_pd/cli/discharge_pd_gene.py",
    "surface_pd/cli/surface_enumeration.py",
    "surface_pd/cli/surface_pd_plot.py",
    "surface_pd/configuration/__init__.py",
    "surface_pd/configuration/phase_diagram.py",
    "surface_pd/core/__init__.py",
    "surface_pd/core/enum.py",
    "surface_pd/core/post_check.py",
    "surface_pd/core/pre_check.py",
    "surface_pd/core/enumeration_slab.py",
    "surface_pd/error/__init__.py",
    "surface_pd/error/error.py",
    "surface_pd/plot/__init__.py",
    "surface_pd/plot/_phase_data_io.py",
    "surface_pd/plot/pd_data.py",
    "surface_pd/plot/plot.py",
    "surface_pd/plot/reference_energies.py",
    "surface_pd/plot/surface_energy.py",
    "surface_pd/schemas/phase-diagram-config-v1.schema.json",
    "surface_pd/util/__init__.py",
    "surface_pd/util/util.py",
)
IMPORT_MODULES = (
    "surface_pd",
    "surface_pd.analysis",
    "surface_pd.cli",
    "surface_pd.configuration",
    "surface_pd.core",
    "surface_pd.error",
    "surface_pd.plot",
    "surface_pd.util",
)


@pytest.fixture(scope="session")
def built_wheel(tmp_path_factory):
    """Build the project wheel into a temporary directory."""
    wheel_dir = tmp_path_factory.mktemp("surface-pd-wheel")
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-deps",
            "--wheel-dir",
            str(wheel_dir),
            ".",
        ],
        cwd=PROJECT_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    wheels = sorted(wheel_dir.glob("surface_pd-*.whl"))
    assert len(wheels) == 1
    return wheels[0]


@pytest.fixture(scope="session")
def installed_wheel_environment(tmp_path_factory, built_wheel):
    """Install the built wheel into an isolated virtual environment."""
    venv_dir = tmp_path_factory.mktemp("surface-pd-install")
    venv.EnvBuilder(with_pip=True, system_site_packages=True).create(venv_dir)

    scripts_dir = venv_dir / ("Scripts" if os.name == "nt" else "bin")
    python = scripts_dir / ("python.exe" if os.name == "nt" else "python")
    result = subprocess.run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--no-deps",
            str(built_wheel),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    return python, scripts_dir


def test_wheel_contains_package_submodules(built_wheel):
    """The built wheel should include all package submodules."""
    with ZipFile(built_wheel) as wheel:
        wheel_names = set(wheel.namelist())

    assert set(PACKAGE_FILES) <= wheel_names
    assert any(
        name.endswith(".dist-info/entry_points.txt")
        for name in wheel_names
    )


def test_wheel_uses_authoritative_package_version(built_wheel):
    """The wheel filename should contain the authoritative release version."""
    assert built_wheel.name.startswith("surface_pd-0.1.0-")


def test_wheel_entry_points_match_console_scripts(built_wheel):
    """Wheel entry-point metadata should expose every console script."""
    with ZipFile(built_wheel) as wheel:
        entry_points_name = next(
            name
            for name in wheel.namelist()
            if name.endswith(".dist-info/entry_points.txt")
        )
        entry_points = wheel.read(entry_points_name).decode()

    for script_name in CONSOLE_SCRIPTS:
        assert script_name in entry_points


def test_installed_wheel_imports_package_modules(installed_wheel_environment):
    """Installed wheels should import the public package submodules."""
    python, _ = installed_wheel_environment
    import_script = "; ".join(f"import {module}" for module in IMPORT_MODULES)
    result = subprocess.run(
        [str(python), "-c", import_script],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_installed_wheel_reports_authoritative_version(
    installed_wheel_environment,
):
    """Installed runtime metadata should match the built wheel version."""
    python, _ = installed_wheel_environment
    result = subprocess.run(
        [
            str(python),
            "-c",
            "import surface_pd; print(surface_pd.__version__)",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "0.1.0"


@pytest.mark.parametrize("script_name", CONSOLE_SCRIPTS)
def test_installed_console_scripts_show_help(
    installed_wheel_environment,
    script_name,
):
    """Installed console scripts should run and expose argparse help."""
    _, scripts_dir = installed_wheel_environment
    script_path = scripts_dir / script_name

    result = subprocess.run(
        [str(script_path), "--help"],
        text=True,
        capture_output=True,
        check=False,
    )

    assert script_path.is_file()
    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout
