"""Read and write self-contained surface phase-diagram data files."""

from io import StringIO
from pathlib import Path

import pandas as pd

from surface_pd.plot.reference_energies import ReferenceEnergies

_METADATA_TO_FIELD = {
    "method": "method",
    "reference_li_ev_per_atom": "li_ev_per_atom",
    "reference_o2_raw_ev_per_molecule": "o2_raw_ev_per_molecule",
    "reference_o2_correction_ev_per_molecule": (
        "o2_correction_ev_per_molecule"
    ),
    "reference_bulk_litmo2_ev_per_formula_unit": (
        "bulk_litmo2_ev_per_formula_unit"
    ),
}


def _parse_metadata(lines: list[str]) -> ReferenceEnergies:
    """Parse and validate a leading phase-data metadata block."""
    metadata = {}
    for line in lines:
        content = line.lstrip()[1:].strip()
        if "=" not in content:
            raise ValueError("Metadata lines must use 'key = value'.")
        key, value = (part.strip() for part in content.split("=", maxsplit=1))
        if key not in _METADATA_TO_FIELD:
            raise ValueError(f"Unknown metadata key {key!r}.")
        if key in metadata:
            raise ValueError(f"Duplicate metadata key {key!r}.")
        metadata[key] = value

    missing = [key for key in _METADATA_TO_FIELD if key not in metadata]
    if missing:
        raise ValueError(
            "Missing required metadata: " + ", ".join(missing) + "."
        )

    values = {"method": metadata["method"]}
    for metadata_key, field_name in _METADATA_TO_FIELD.items():
        if metadata_key == "method":
            continue
        try:
            values[field_name] = float(metadata[metadata_key])
        except ValueError as error:
            raise ValueError(
                f"{metadata_key} must be a finite number."
            ) from error
    return ReferenceEnergies(**values)


def _read_phase_diagram_file(
    path: str | Path,
) -> tuple[pd.DataFrame, ReferenceEnergies]:
    """Read a phase table and its required reference-energy metadata."""
    lines = Path(path).read_text().splitlines(keepends=True)
    metadata_lines = []
    table_start = None
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if table_start is not None:
                raise ValueError("Metadata must precede the tabular header.")
            metadata_lines.append(line)
            continue
        table_start = index
        break

    references = _parse_metadata(metadata_lines)
    if table_start is None:
        raise ValueError("Phase-data file must contain a tabular header.")
    table_lines = lines[table_start:]
    if any(line.strip().startswith("#") for line in table_lines[1:]):
        raise ValueError("Metadata must precede the tabular header.")

    try:
        dataframe = pd.read_csv(
            StringIO("".join(table_lines)),
            sep=r"\s+",
            index_col=0,
        )
    except (pd.errors.ParserError, pd.errors.EmptyDataError) as error:
        raise ValueError("Phase-data table is malformed.") from error
    if dataframe.empty:
        raise ValueError(
            "Phase-data table must contain at least one phase row."
        )
    return dataframe, references


def _write_phase_diagram_file(
    path: str | Path,
    dataframe: pd.DataFrame,
    references: ReferenceEnergies,
) -> None:
    """Write phase data with canonical, complete reference metadata."""
    metadata = (
        f"# method = {references.method}\n"
        f"# reference_li_ev_per_atom = {references.li_ev_per_atom}\n"
        "# reference_o2_raw_ev_per_molecule = "
        f"{references.o2_raw_ev_per_molecule}\n"
        "# reference_o2_correction_ev_per_molecule = "
        f"{references.o2_correction_ev_per_molecule}\n"
        "# reference_bulk_litmo2_ev_per_formula_unit = "
        f"{references.bulk_litmo2_ev_per_formula_unit}\n"
    )
    table = dataframe.to_csv(sep="\t")
    Path(path).write_text(metadata + table)


def _require_matching_reference_energies(
    first: ReferenceEnergies,
    second: ReferenceEnergies,
) -> None:
    """Reject reference metadata that differs between aligned datasets."""
    if first == second:
        return
    field_names = (
        "method",
        "li_ev_per_atom",
        "o2_raw_ev_per_molecule",
        "o2_correction_ev_per_molecule",
        "bulk_litmo2_ev_per_formula_unit",
    )
    differences = [
        field_name
        for field_name in field_names
        if getattr(first, field_name) != getattr(second, field_name)
    ]
    raise ValueError(
        "Input files use incompatible reference energies; differing fields: "
        + ", ".join(differences)
        + "."
    )
