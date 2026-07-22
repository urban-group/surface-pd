"""Strict JSON configuration for generalized phase diagrams."""

import json
from collections.abc import Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

from surface_pd.thermodynamics import (
    AlignedPhaseDataset,
    ConstantChemicalPotential,
    DatasetAlignment,
    DiagramAxis,
    DirectChemicalPotential,
    FixedPressureOxygenChemicalPotential,
    GrandPotentialModel,
    IntercalationChemicalPotential,
    PhaseDataset,
    PhaseDiagramSpecification,
    ReferencePhase,
)
from surface_pd.thermodynamics._validation import (
    validate_identifier,
    validate_positive_integer,
    validate_provenance,
    validate_variable_name,
)

from ._table_loading import load_phase_dataset

_SCHEMA_VERSION = 1
_MODEL_NAMES = {
    "constant",
    "direct",
    "intercalation_voltage",
    "fixed_pressure_oxygen",
}


def _require_mapping(value: object, context: str) -> dict[str, Any]:
    """Return a plain mapping or raise a contextual type error."""
    if not isinstance(value, Mapping):
        raise TypeError(f"{context} must be a JSON object")
    if any(not isinstance(key, str) for key in value):
        raise TypeError(f"{context} field names must be strings")
    return dict(value)


def _require_keys(
    value: object,
    required: set[str],
    context: str,
    optional: set[str] | None = None,
) -> dict[str, Any]:
    """Return a mapping with exactly its allowed versioned fields."""
    mapping = _require_mapping(value, context)
    optional = optional or set()
    missing = sorted(required - set(mapping))
    unknown = sorted(set(mapping) - required - optional)
    if missing:
        raise ValueError(f"{context} has missing fields: {', '.join(missing)}")
    if unknown:
        raise ValueError(f"{context} has unknown fields: {', '.join(unknown)}")
    return mapping


def _require_sequence(value: object, context: str) -> list[Any]:
    """Return a JSON array as a list."""
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{context} must be a JSON array")
    return list(value)


def _require_column_name(value: object, context: str) -> str:
    """Return a nonempty single-line table-column name."""
    return validate_provenance(value, context)


def _parse_model(value: object, component: str) -> object:
    """Construct one reviewed built-in chemical-potential model."""
    context = f"chemical-potential model for {component!r}"
    mapping = _require_mapping(value, context)
    model_name = mapping.get("model")
    if model_name not in _MODEL_NAMES:
        raise ValueError(f"unknown chemical-potential model {model_name!r}")

    if model_name == "constant":
        mapping = _require_keys(
            mapping,
            {"model", "value_ev_per_component"},
            "constant model",
        )
        return ConstantChemicalPotential(mapping["value_ev_per_component"])
    if model_name == "direct":
        mapping = _require_keys(
            mapping,
            {
                "model",
                "state_variable",
                "reference_energy_ev_per_component",
            },
            "direct model",
        )
        return DirectChemicalPotential(
            mapping["state_variable"],
            mapping["reference_energy_ev_per_component"],
        )
    if model_name == "intercalation_voltage":
        mapping = _require_keys(
            mapping,
            {
                "model",
                "reference_energy_ev_per_component",
                "electrons_per_component",
                "voltage_reference",
            },
            "intercalation_voltage model",
        )
        return IntercalationChemicalPotential(
            mapping["reference_energy_ev_per_component"],
            mapping["electrons_per_component"],
            mapping["voltage_reference"],
        )

    mapping = _require_keys(
        mapping,
        {
            "model",
            "raw_o2_energy_ev_per_molecule",
            "correction_ev_per_molecule",
        },
        "fixed_pressure_oxygen model",
    )
    return FixedPressureOxygenChemicalPotential(
        mapping["raw_o2_energy_ev_per_molecule"],
        mapping["correction_ev_per_molecule"],
    )


def _parse_axis_metadata(
    value: object, context: str
) -> tuple[str, str, str]:
    """Return one configured state-variable identity and display metadata."""
    mapping = _require_keys(
        value,
        {"state_variable", "label", "unit"},
        context,
    )
    return (
        validate_variable_name(mapping["state_variable"]),
        validate_provenance(mapping["label"], f"{context} label"),
        validate_provenance(mapping["unit"], f"{context} unit"),
    )


def _resolved_columns(
    dataset: Mapping[str, Any], components: tuple[str, ...]
) -> dict[str, Any]:
    """Return canonical table columns with validated optional overrides."""
    columns: dict[str, Any] = {
        "phase_id": "phase_id",
        "composition": {component: component for component in components},
        "dft_energy_ev": "dft_energy_ev",
        "surface_area_angstrom2": "surface_area_angstrom2",
    }
    overrides = _require_keys(
        dataset.get("column_overrides", {}),
        set(),
        "column_overrides",
        {
            "phase_id",
            "composition",
            "dft_energy_ev",
            "surface_area_angstrom2",
        },
    )
    for field_name in (
        "phase_id",
        "dft_energy_ev",
        "surface_area_angstrom2",
    ):
        if field_name in overrides:
            columns[field_name] = _require_column_name(
                overrides[field_name], f"{field_name} column override"
            )
    if "composition" in overrides:
        composition = _require_mapping(
            overrides["composition"], "composition column overrides"
        )
        unknown = sorted(set(composition) - set(components))
        if unknown:
            raise ValueError(
                "composition column overrides contain unknown components: "
                + ", ".join(unknown)
            )
        for component, column in composition.items():
            columns["composition"][component] = _require_column_name(
                column, f"composition column override {component}"
            )
    mapped_columns = [
        columns["phase_id"],
        columns["dft_energy_ev"],
        columns["surface_area_angstrom2"],
        *columns["composition"].values(),
    ]
    if len(mapped_columns) != len(set(mapped_columns)):
        raise ValueError("resolved source columns must be unique")
    return columns


class PhaseDiagramConfiguration:
    """Own a validated version-1 generalized phase-diagram configuration.

    Parameters
    ----------
    data : mapping
        Complete JSON-compatible version-1 configuration.
    source_path : path-like or None, optional
        Source JSON path retained for future relative dataset-path resolution.

    Notes
    -----
    JSON selects only reviewed package model names. It cannot import arbitrary
    Python objects or serialize custom callables.
    """

    __slots__ = (
        "_calculation_method",
        "_components",
        "_data",
        "_model",
        "_source_path",
        "_x_axis_metadata",
        "_y_axis_metadata",
    )

    def __init__(
        self,
        data: Mapping[str, Any],
        *,
        source_path: str | Path | None = None,
    ) -> None:
        try:
            canonical = json.loads(json.dumps(data, allow_nan=False))
        except (TypeError, ValueError) as error:
            raise TypeError(
                "configuration data must contain only JSON-compatible values"
            ) from error
        self._validate_and_construct(canonical)
        self._data = canonical
        self._source_path = (
            None if source_path is None else Path(source_path).resolve()
        )

    def _validate_and_construct(self, data: object) -> None:
        """Validate the full schema and construct its domain objects."""
        root = _require_keys(
            data,
            {
                "schema_version",
                "calculation_method",
                "components",
                "independent_chemical_potentials",
                "reference_phases",
                "diagram",
                "datasets",
                "alignments",
            },
            "configuration",
        )
        version = root["schema_version"]
        if isinstance(version, bool) or version != _SCHEMA_VERSION:
            raise ValueError(f"unsupported schema_version {version!r}")
        method = validate_provenance(
            root["calculation_method"], "calculation_method"
        )
        components = tuple(
            validate_variable_name(component)
            for component in _require_sequence(
                root["components"], "components"
            )
        )
        if not components or len(set(components)) != len(components):
            raise ValueError("components must be nonempty and unique")

        model_entries = _require_mapping(
            root["independent_chemical_potentials"],
            "independent_chemical_potentials",
        )
        unknown_components = set(model_entries) - set(components)
        if unknown_components:
            raise ValueError(
                "independent chemical potentials contain unknown components: "
                + ", ".join(sorted(unknown_components))
            )
        independent_models = {
            component: _parse_model(model_entries[component], component)
            for component in components
            if component in model_entries
        }

        references = []
        for index, item in enumerate(
            _require_sequence(root["reference_phases"], "reference_phases")
        ):
            reference = _require_keys(
                item,
                {
                    "reference_id",
                    "composition",
                    "energy_ev_per_formula_unit",
                },
                f"reference_phases[{index}]",
            )
            references.append(
                ReferencePhase(
                    reference["reference_id"],
                    reference["composition"],
                    reference["energy_ev_per_formula_unit"],
                    method,
                )
            )
        model = GrandPotentialModel(
            components, independent_models, tuple(references)
        )

        diagram_data = _require_keys(
            root["diagram"],
            {"x_axis", "y_axis"},
            "diagram",
        )
        x_axis_metadata = _parse_axis_metadata(
            diagram_data["x_axis"], "x_axis"
        )
        y_axis_metadata = _parse_axis_metadata(
            diagram_data["y_axis"], "y_axis"
        )
        axis_variables = {x_axis_metadata[0], y_axis_metadata[0]}
        if len(axis_variables) != 2:
            raise ValueError(
                "x_axis and y_axis state variables must be distinct"
            )
        unknown_axis_variables = axis_variables - set(
            model.required_state_variables
        )
        if unknown_axis_variables:
            raise ValueError(
                "diagram axes contain unused state variables: "
                + ", ".join(sorted(unknown_axis_variables))
            )

        dataset_ids = self._validate_datasets(root["datasets"], components)
        reference_ids = {reference.reference_id for reference in references}
        self._validate_alignments(
            root["alignments"], dataset_ids, reference_ids
        )
        self._calculation_method = method
        self._components = components
        self._model = model
        self._x_axis_metadata = x_axis_metadata
        self._y_axis_metadata = y_axis_metadata

    @staticmethod
    def _validate_datasets(
        value: object, components: tuple[str, ...]
    ) -> set[str]:
        """Validate explicit phase-table definitions and return their IDs."""
        items = _require_sequence(value, "datasets")
        if not items:
            raise ValueError("datasets must contain at least one definition")
        dataset_ids = set()
        for index, item in enumerate(items):
            dataset = _require_keys(
                item,
                {"dataset_id", "path", "number_of_surfaces"},
                f"datasets[{index}]",
                {"column_overrides"},
            )
            dataset_id = validate_identifier(
                dataset["dataset_id"], f"datasets[{index}].dataset_id"
            )
            if dataset_id in dataset_ids:
                raise ValueError("dataset_id values must be unique")
            dataset_ids.add(dataset_id)
            _require_column_name(dataset["path"], f"datasets[{index}].path")
            validate_positive_integer(
                dataset["number_of_surfaces"], "number_of_surfaces"
            )
            _resolved_columns(dataset, components)
        return dataset_ids

    @staticmethod
    def _validate_alignments(
        value: object,
        dataset_ids: set[str],
        reference_ids: set[str],
    ) -> None:
        """Validate direct-to-reference alignment identities."""
        targets = set()
        for index, item in enumerate(
            _require_sequence(value, "alignments")
        ):
            alignment = _require_keys(
                item,
                {
                    "reference_dataset_id",
                    "target_dataset_id",
                    "reference_anchor_phase_id",
                    "target_anchor_phase_id",
                    "bulk_reference_id",
                },
                f"alignments[{index}]",
            )
            reference_dataset_id = alignment["reference_dataset_id"]
            target_id = alignment["target_dataset_id"]
            if (
                reference_dataset_id not in dataset_ids
                or target_id not in dataset_ids
            ):
                raise ValueError("alignment contains an unknown dataset_id")
            if reference_dataset_id == target_id:
                raise ValueError(
                    "alignment reference and target must be distinct"
                )
            if target_id in targets:
                raise ValueError("each alignment target must be unique")
            targets.add(target_id)
            validate_identifier(
                alignment["reference_anchor_phase_id"],
                "reference_anchor_phase_id",
            )
            validate_identifier(
                alignment["target_anchor_phase_id"],
                "target_anchor_phase_id",
            )
            if alignment["bulk_reference_id"] not in reference_ids:
                raise ValueError(
                    "alignment contains an unknown bulk_reference_id"
                )
        reference_dataset_ids = {
            alignment["reference_dataset_id"]
            for alignment in _require_sequence(value, "alignments")
        }
        if reference_dataset_ids & targets:
            raise ValueError(
                "an alignment reference cannot also be an alignment target"
            )

    @property
    def schema_version(self) -> int:
        """Return the supported configuration schema version."""
        return _SCHEMA_VERSION

    @property
    def calculation_method(self) -> str:
        """Return the uninterpreted calculation-method provenance."""
        return self._calculation_method

    @property
    def components(self) -> tuple[str, ...]:
        """Return the ordered thermodynamic component basis."""
        return self._components

    @property
    def model(self) -> GrandPotentialModel:
        """Return the validated generalized grand-potential model."""
        return self._model

    def create_diagram_specification(
        self,
        *,
        x_values: object,
        y_values: object,
        fixed_conditions: Mapping[str, object] | None = None,
    ) -> PhaseDiagramSpecification:
        """Construct one explicit numerical evaluation domain.

        Parameters
        ----------
        x_values, y_values : array-like
            Finite, strictly increasing coordinates supplied for this
            evaluation. Each axis requires at least two values.
        fixed_conditions : mapping or None, optional
            Finite scalar values for required model state variables that are
            not assigned to either configured axis.

        Returns
        -------
        PhaseDiagramSpecification
            Validated numerical specification using the configured state-
            variable identities and default display metadata.
        """
        fixed = {} if fixed_conditions is None else fixed_conditions
        specification = PhaseDiagramSpecification(
            DiagramAxis(
                self._x_axis_metadata[0],
                x_values,
                self._x_axis_metadata[1],
                self._x_axis_metadata[2],
            ),
            DiagramAxis(
                self._y_axis_metadata[0],
                y_values,
                self._y_axis_metadata[1],
                self._y_axis_metadata[2],
            ),
            fixed,
        )
        supplied = {
            specification.x_axis.state_variable,
            specification.y_axis.state_variable,
            *specification.fixed_conditions,
        }
        required = set(self.model.required_state_variables)
        if supplied != required:
            missing = sorted(required - supplied)
            unused = sorted(supplied - required)
            details = []
            if missing:
                details.append("missing: " + ", ".join(missing))
            if unused:
                details.append("unused: " + ", ".join(unused))
            raise ValueError(
                "diagram state variables do not match model requirements; "
                + "; ".join(details)
            )
        return specification

    @property
    def source_path(self) -> Path | None:
        """Return the absolute source JSON path, if read from a file."""
        return self._source_path

    def to_dict(self) -> dict[str, Any]:
        """Return a deep mutable copy of the canonical JSON data."""
        return deepcopy(self._data)

    def load_datasets(
        self,
    ) -> tuple[PhaseDataset | AlignedPhaseDataset, ...]:
        """Load configured phase tables and apply direct alignments.

        Returns
        -------
        tuple of PhaseDataset or AlignedPhaseDataset
            Immutable datasets in configuration order. Declared alignment
            targets are returned as non-mutating aligned views.

        Raises
        ------
        ValueError
            If a relative table path has no source JSON directory, a table is
            malformed, mapped data are invalid, or alignment is incompatible.
        """
        datasets = {}
        for definition in self._data["datasets"]:
            path = Path(validate_provenance(definition["path"], "path"))
            if not path.is_absolute():
                if self._source_path is None:
                    raise ValueError(
                        "relative dataset paths require a source JSON path"
                    )
                path = self._source_path.parent / path
            dataset = load_phase_dataset(
                definition,
                path=path,
                components=self._components,
                calculation_method=self._calculation_method,
            )
            datasets[dataset.dataset_id] = dataset

        references = {
            reference.reference_id: reference
            for reference in self._model.reference_phases
        }
        aligned = {}
        for definition in self._data["alignments"]:
            reference_dataset_id = validate_identifier(
                definition["reference_dataset_id"], "reference_dataset_id"
            )
            target_id = validate_identifier(
                definition["target_dataset_id"], "target_dataset_id"
            )
            bulk_reference_id = validate_identifier(
                definition["bulk_reference_id"], "bulk_reference_id"
            )
            alignment = DatasetAlignment(
                datasets[reference_dataset_id],
                datasets[target_id],
                definition["reference_anchor_phase_id"],
                definition["target_anchor_phase_id"],
                references[bulk_reference_id],
            )
            aligned[target_id] = alignment.create_aligned_dataset()

        return tuple(
            aligned.get(dataset_id, dataset)
            for dataset_id, dataset in datasets.items()
        )

    @classmethod
    def read_json(cls, path: str | Path) -> "PhaseDiagramConfiguration":
        """Read and validate one UTF-8 JSON configuration file.

        Parameters
        ----------
        path : path-like
            Configuration file to read.
        """
        source = Path(path)
        try:
            data = json.loads(source.read_text())
        except json.JSONDecodeError as error:
            raise ValueError(
                f"{source} must contain valid JSON: {error.msg}"
            ) from error
        return cls(data, source_path=source)

    def write_json(self, path: str | Path) -> None:
        """Write canonical version-1 JSON with stable formatting.

        Parameters
        ----------
        path : path-like
            Destination file, which is replaced if it already exists.
        """
        destination = Path(path)
        destination.write_text(json.dumps(self._data, indent=2) + "\n")
