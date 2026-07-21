"""Strict loading of canonical version-1 phase tables."""

from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from surface_pd.thermodynamics import Phase, PhaseDataset
from surface_pd.thermodynamics._validation import validate_identifier


def _finite_float(value: str, context: str) -> float:
    """Parse one finite floating-point table field."""
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{context} must be a finite number") from error
    if not np.isfinite(result):
        raise ValueError(f"{context} must be a finite number")
    return result


def _integer(value: str, context: str, *, positive: bool) -> int:
    """Parse an exact integer, optionally requiring strict positivity."""
    qualifier = "positive " if positive else "nonnegative "
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"{context} must be a {qualifier}integer") from error
    if not np.isfinite(number):
        raise ValueError(f"{context} must be a {qualifier}integer")
    if not number.is_integer():
        raise ValueError(f"{context} must be a {qualifier}integer")
    result = int(number)
    if (positive and result <= 0) or (not positive and result < 0):
        raise ValueError(f"{context} must be a {qualifier}integer")
    return result


def _read_rows(
    path: Path, dataset_id: str
) -> tuple[list[str], list[tuple[int, list[str]]]]:
    """Read a strict whitespace table while retaining physical line numbers."""
    try:
        lines = path.read_text().splitlines()
    except OSError as error:
        raise ValueError(
            f"dataset {dataset_id!r} could not read table {path}: {error}"
        ) from error

    content = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            raise ValueError(
                f"dataset {dataset_id!r}: comments are not supported in "
                "version-1 phase tables"
            )
        content.append((line_number, stripped.split()))
    if not content:
        raise ValueError(f"dataset {dataset_id!r} table is empty")

    _, header = content[0]
    if len(set(header)) != len(header):
        raise ValueError(
            f"dataset {dataset_id!r} table header columns must be unique"
        )
    rows = content[1:]
    if not rows:
        raise ValueError(
            f"dataset {dataset_id!r} table must contain at least one phase row"
        )
    for line_number, fields in rows:
        if len(fields) != len(header):
            raise ValueError(
                f"dataset {dataset_id!r} row {line_number} has "
                f"{len(fields)} fields; expected {len(header)} fields"
            )
    return header, rows


def load_phase_dataset(
    definition: Mapping[str, object],
    *,
    path: Path,
    components: Sequence[str],
    calculation_method: str,
) -> PhaseDataset:
    """Load one configured table as an immutable phase dataset."""
    dataset_id = validate_identifier(definition["dataset_id"], "dataset_id")
    columns = {
        "phase_id": "phase_id",
        "composition": {component: component for component in components},
        "dft_energy_ev": "dft_energy_ev",
        "surface_area_angstrom2": "surface_area_angstrom2",
    }
    overrides = definition.get("column_overrides", {})
    for field_name in (
        "phase_id",
        "dft_energy_ev",
        "surface_area_angstrom2",
    ):
        if field_name in overrides:
            columns[field_name] = overrides[field_name]
    columns["composition"].update(overrides.get("composition", {}))
    header, rows = _read_rows(path, dataset_id)

    mapped_columns = {
        columns["phase_id"],
        columns["dft_energy_ev"],
        *columns["composition"].values(),
    }
    mapped_columns.add(columns["surface_area_angstrom2"])
    missing = sorted(mapped_columns - set(header))
    if missing:
        raise ValueError(
            f"dataset {dataset_id!r} has missing required columns: "
            + ", ".join(missing)
        )

    phases = []
    phase_ids = set()
    for line_number, fields in rows:
        row = dict(zip(header, fields, strict=True))
        phase_id = row[columns["phase_id"]]
        if phase_id in phase_ids:
            raise ValueError(f"duplicate phase_id {phase_id!r}")
        phase_ids.add(phase_id)
        composition = {
            component: _integer(
                row[columns["composition"][component]],
                f"dataset {dataset_id!r} row {line_number} column "
                f"{columns['composition'][component]!r}",
                positive=False,
            )
            for component in components
        }
        energy_column = columns["dft_energy_ev"]
        energy = _finite_float(
            row[energy_column],
            f"dataset {dataset_id!r} row {line_number} column "
            f"{energy_column!r}",
        )
        area_column = columns["surface_area_angstrom2"]
        area = _finite_float(
            row[area_column],
            f"dataset {dataset_id!r} row {line_number} column "
            f"{area_column!r}",
        )
        number_of_surfaces = definition["number_of_surfaces"]
        try:
            phases.append(
                    Phase(
                        phase_id,
                        composition,
                        energy,
                        area,
                        number_of_surfaces,
                    )
            )
        except (TypeError, ValueError) as error:
            raise type(error)(
                f"dataset {dataset_id!r} row {line_number}: {error}"
            ) from error

    return PhaseDataset(
        dataset_id,
        components,
        phases,
        calculation_method,
    )
