"""Contracts preventing return of the chemistry-specific legacy workflow."""

import importlib.util
import tomllib
from pathlib import Path

from surface_pd.thermodynamics import (
    ConstantChemicalPotential,
    GrandPotentialModel,
    Phase,
    PhaseDataset,
    ThermodynamicState,
)

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_REMOVED_MODULES = (
    "surface_pd.cli.discharge_pd_gene",
    "surface_pd.plot._phase_data_io",
    "surface_pd.plot.pd_data",
    "surface_pd.plot.plot",
    "surface_pd.plot.reference_energies",
    "surface_pd.plot.surface_energy",
)


def test_chemistry_specific_legacy_modules_are_absent():
    """Removed pre-1.0 implementations should not remain importable."""
    for module_name in _REMOVED_MODULES:
        assert importlib.util.find_spec(module_name) is None


def test_discharge_generator_console_script_is_absent():
    """Candidate filtering should not be exposed as package policy."""
    metadata = tomllib.loads((_PROJECT_ROOT / "pyproject.toml").read_text())

    assert "generate-discharge-pd" not in metadata["project"]["scripts"]


def test_zero_dft_energy_is_an_ordinary_thermodynamic_value():
    """Zero must be evaluated numerically rather than treated as exclusion."""
    phase = Phase("zero", {"A": 1}, 0.0, 2.0, 1)
    dataset = PhaseDataset("data", ("A",), (phase,), "test method")
    model = GrandPotentialModel(
        ("A",), {"A": ConstantChemicalPotential(1.5)}, ()
    )

    result = model.evaluate(dataset, ThermodynamicState({}))

    assert result.total_grand_potential_ev.tolist() == [-1.5]
    assert result.surface_grand_potential_ev_per_angstrom2.tolist() == [-0.75]
