"""Rendering for generalized numerical phase diagrams."""

from dataclasses import dataclass
from numbers import Real
from typing import Literal

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.colorbar import Colorbar
from matplotlib.colors import BoundaryNorm, is_color_like
from matplotlib.figure import Figure
from matplotlib.typing import ColorType

from surface_pd.thermodynamics import PhaseDiagramResult
from surface_pd.thermodynamics._validation import (
    validate_provenance,
    validate_variable_name,
)


@dataclass(frozen=True)
class CompositionColoring:
    """Define an explicit phase-composition color scale.

    Parameters
    ----------
    component : str
        Composition component used as the numerator.
    normalization : {"atomic_fraction", "component_ratio"}
        ``"atomic_fraction"`` calculates the component count divided by the
        sum of all component counts. ``"component_ratio"`` divides by the
        explicitly selected ``reference_component`` count.
    reference_component : str or None
        Denominator component for ``"component_ratio"``; it must be ``None``
        for ``"atomic_fraction"``.
    label : str
        Nonempty, single-line colorbar label without units.
    unit : str
        Nonempty, single-line display unit.
    """

    component: str
    normalization: Literal["atomic_fraction", "component_ratio"]
    reference_component: str | None
    label: str
    unit: str

    def __post_init__(self) -> None:
        """Validate the complete composition-normalization definition."""
        component = validate_variable_name(self.component)
        if self.normalization not in {
            "atomic_fraction",
            "component_ratio",
        }:
            raise ValueError(
                "normalization must be 'atomic_fraction' or "
                "'component_ratio'"
            )
        if self.normalization == "atomic_fraction":
            if self.reference_component is not None:
                raise ValueError(
                    "reference_component must be None for atomic_fraction"
                )
            reference_component = None
        else:
            if self.reference_component is None:
                raise ValueError(
                    "reference_component is required for component_ratio"
                )
            reference_component = validate_variable_name(
                self.reference_component
            )
            if reference_component == component:
                raise ValueError(
                    "component and reference_component must be distinct"
                )

        object.__setattr__(self, "component", component)
        object.__setattr__(
            self, "reference_component", reference_component
        )
        object.__setattr__(
            self, "label", validate_provenance(self.label, "label")
        )
        object.__setattr__(
            self, "unit", validate_provenance(self.unit, "unit")
        )

    def phase_values(self, result: PhaseDiagramResult) -> np.ndarray:
        """Return one composition value per phase in result order.

        Parameters
        ----------
        result : PhaseDiagramResult
            Numerical diagram whose candidate phase compositions are used.

        Returns
        -------
        numpy.ndarray
            Read-only values in the same order as ``result.phase_ids``.
        """
        if not isinstance(result, PhaseDiagramResult):
            raise TypeError("result must be a PhaseDiagramResult")

        values = []
        for dataset in result.datasets:
            if self.component not in dataset.components:
                raise ValueError(
                    f"composition component {self.component!r} is not present "
                    f"in dataset {dataset.dataset_id!r}"
                )
            if (
                self.reference_component is not None
                and self.reference_component not in dataset.components
            ):
                raise ValueError(
                    "reference component "
                    f"{self.reference_component!r} is not present in dataset "
                    f"{dataset.dataset_id!r}"
                )
            for phase in dataset.phases:
                numerator = phase.composition[self.component]
                if self.normalization == "atomic_fraction":
                    denominator = sum(phase.composition.values())
                else:
                    denominator = phase.composition[
                        self.reference_component
                    ]
                    if denominator == 0:
                        raise ValueError(
                            "reference component "
                            f"{self.reference_component!r} has zero count in "
                            "phase "
                            f"{dataset.qualified_phase_id(phase.phase_id)!r}"
                        )
                values.append(numerator / denominator)

        result_values = np.array(values, dtype=float)
        result_values.setflags(write=False)
        return result_values


def _axis_label(label: str, unit: str) -> str:
    """Return one explicit Matplotlib axis or colorbar label."""
    return f"{label} ({unit})"


def _cell_edges(coordinates: np.ndarray) -> np.ndarray:
    """Return nearest-shading cell edges for increasing center coordinates."""
    midpoints = 0.5 * (coordinates[:-1] + coordinates[1:])
    first = coordinates[0] - 0.5 * (coordinates[1] - coordinates[0])
    last = coordinates[-1] + 0.5 * (coordinates[-1] - coordinates[-2])
    return np.concatenate(([first], midpoints, [last]))


def _phase_boundary_segments(result: PhaseDiagramResult) -> list[np.ndarray]:
    """Return cell-edge segments separating representative stable phases."""
    representatives = result.representative_phase_indices
    x_values = result.specification.x_axis.values
    y_values = result.specification.y_axis.values
    x_edges = _cell_edges(x_values)
    y_edges = _cell_edges(y_values)
    segments = []

    for y_index, x_index in np.argwhere(
        representatives[:, 1:] != representatives[:, :-1]
    ):
        x_boundary = x_edges[x_index + 1]
        segments.append(
            np.array(
                [
                    [x_boundary, y_edges[y_index]],
                    [x_boundary, y_edges[y_index + 1]],
                ]
            )
        )
    for y_index, x_index in np.argwhere(
        representatives[1:, :] != representatives[:-1, :]
    ):
        y_boundary = y_edges[y_index + 1]
        segments.append(
            np.array(
                [
                    [x_edges[x_index], y_boundary],
                    [x_edges[x_index + 1], y_boundary],
                ]
            )
        )
    return segments


def plot_phase_diagram(
    result: PhaseDiagramResult,
    *,
    coloring: CompositionColoring | None = None,
    ax: Axes | None = None,
    cmap: str | None = None,
    invert_x_axis: bool = False,
    invert_y_axis: bool = False,
    boundary_color: ColorType | None = "black",
    boundary_linewidth: float = 1.0,
) -> tuple[Figure, Axes, Colorbar]:
    """Render an already evaluated generalized phase diagram.

    Parameters
    ----------
    result : PhaseDiagramResult
        Complete numerical result. Rendering does not reevaluate it.
    coloring : CompositionColoring or None, optional
        Explicit continuous composition scale. ``None`` uses discrete
        qualified phase identities.
    ax : matplotlib.axes.Axes or None, optional
        Existing axes to draw into. New figure and axes are created when
        omitted.
    cmap : str or None, optional
        Matplotlib colormap name. Defaults to ``"tab20"`` for phase identity
        and ``"viridis"`` for composition values.
    invert_x_axis : bool, optional
        If true, present the x-axis in decreasing screen order.
    invert_y_axis : bool, optional
        If true, present the y-axis in decreasing screen order.
    boundary_color : color-like or None, optional
        Matplotlib color for boundaries between representative stable phases.
        Defaults to black. ``None`` disables the boundary overlay.
    boundary_linewidth : float, optional
        Positive finite phase-boundary width in points. This option is ignored
        when ``boundary_color`` is ``None``.

    Returns
    -------
    tuple
        Matplotlib ``(figure, axes, colorbar)`` objects. The function neither
        displays nor saves the figure.
    """
    if not isinstance(result, PhaseDiagramResult):
        raise TypeError("result must be a PhaseDiagramResult")
    if coloring is not None and not isinstance(coloring, CompositionColoring):
        raise TypeError("coloring must be a CompositionColoring or None")
    if ax is not None and not isinstance(ax, Axes):
        raise TypeError("ax must be a matplotlib.axes.Axes or None")
    if not isinstance(invert_x_axis, bool) or not isinstance(
        invert_y_axis, bool
    ):
        raise TypeError("axis inversion options must be boolean")
    if boundary_color is not None and not is_color_like(boundary_color):
        raise ValueError("boundary_color must be a Matplotlib color or None")
    if boundary_color is not None and (
        isinstance(boundary_linewidth, bool)
        or not isinstance(boundary_linewidth, Real)
        or not np.isfinite(boundary_linewidth)
        or boundary_linewidth <= 0
    ):
        raise ValueError(
            "boundary_linewidth must be a positive finite real number"
        )

    if ax is None:
        figure, axes = plt.subplots()
    else:
        axes = ax
        figure = axes.figure

    representatives = result.representative_phase_indices
    if coloring is None:
        phase_count = len(result.phase_ids)
        color_values = representatives
        color_map = matplotlib.colormaps.get_cmap(
            cmap or "tab20"
        ).resampled(phase_count)
        norm = BoundaryNorm(
            np.arange(phase_count + 1, dtype=float) - 0.5,
            color_map.N,
        )
        mesh = axes.pcolormesh(
            result.x_mesh,
            result.y_mesh,
            color_values,
            shading="nearest",
            cmap=color_map,
            norm=norm,
        )
        present_indices = np.unique(representatives)
        colorbar = figure.colorbar(mesh, ax=axes, ticks=present_indices)
        colorbar.ax.set_yticklabels(
            [result.phase_ids[index] for index in present_indices]
        )
        colorbar.set_label("Stable phase")
    else:
        phase_values = coloring.phase_values(result)
        color_values = phase_values[representatives]
        mesh = axes.pcolormesh(
            result.x_mesh,
            result.y_mesh,
            color_values,
            shading="nearest",
            cmap=cmap or "viridis",
        )
        colorbar = figure.colorbar(mesh, ax=axes)
        colorbar.set_label(_axis_label(coloring.label, coloring.unit))

    if boundary_color is not None:
        boundary_segments = _phase_boundary_segments(result)
        if boundary_segments:
            axes.add_collection(
                LineCollection(
                    boundary_segments,
                    colors=boundary_color,
                    linewidths=float(boundary_linewidth),
                    zorder=mesh.get_zorder() + 1,
                )
            )

    x_axis = result.specification.x_axis
    y_axis = result.specification.y_axis
    axes.set_xlabel(_axis_label(x_axis.label, x_axis.unit))
    axes.set_ylabel(_axis_label(y_axis.label, y_axis.unit))
    if invert_x_axis and not axes.xaxis_inverted():
        axes.invert_xaxis()
    if invert_y_axis and not axes.yaxis_inverted():
        axes.invert_yaxis()
    return figure, axes, colorbar
