#!/usr/bin/env python3

"""Plot a generalized phase diagram from versioned JSON configuration."""

import argparse
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.figure import Figure
from matplotlib.legend import Legend

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.plot import CompositionColoring, plot_phase_diagram
from surface_pd.thermodynamics import PhaseDiagramResult


def _fixed_condition(value: str) -> tuple[str, float]:
    """Parse one ``STATE_VARIABLE=VALUE`` command-line condition."""
    variable, separator, raw_value = value.partition("=")
    if not separator or not variable or not raw_value:
        raise argparse.ArgumentTypeError(
            "condition must have the form STATE_VARIABLE=VALUE"
        )
    try:
        return variable, float(raw_value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"fixed condition {variable!r} must be numerical"
        ) from error


def _prepare_phase_diagram(
    configuration: PhaseDiagramConfiguration,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    mesh_points: int = 201,
    fixed_conditions: dict[str, float] | None = None,
    color_component: str | None = None,
    phase_identity_colors: bool = False,
) -> tuple[PhaseDiagramResult, Figure, Axes, Colorbar | Legend]:
    """Evaluate and render without displaying or saving the figure."""
    if not isinstance(configuration, PhaseDiagramConfiguration):
        raise TypeError(
            "configuration must be a PhaseDiagramConfiguration"
        )
    datasets = configuration.load_datasets()
    if isinstance(mesh_points, bool) or not isinstance(mesh_points, int):
        raise TypeError("mesh_points must be an integer")
    if mesh_points < 2:
        raise ValueError("mesh_points must be at least 2")
    specification = configuration.create_diagram_specification(
        x_values=np.linspace(*x_range, mesh_points),
        y_values=np.linspace(*y_range, mesh_points),
        fixed_conditions=fixed_conditions,
    )
    result = specification.evaluate(
        configuration.model,
        datasets,
    )
    coloring = None
    if phase_identity_colors:
        coloring = "phase_identity"
    elif color_component is not None:
        if color_component not in result.independent_components:
            raise ValueError(
                "color component must be an independent component: "
                + ", ".join(result.independent_components)
            )
        coloring = CompositionColoring.atomic_fraction(color_component)
    figure, axes, color_guide = plot_phase_diagram(
        result,
        coloring=coloring,
    )
    return result, figure, axes, color_guide


def _parser() -> argparse.ArgumentParser:
    """Return the generalized plotting argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate and plot a generalized two-dimensional phase diagram "
            "from versioned JSON configuration."
        )
    )
    parser.add_argument(
        "configuration",
        metavar="CONFIG",
        type=Path,
        help="path to a versioned phase-diagram JSON configuration",
    )
    parser.add_argument(
        "--x-range",
        nargs=2,
        type=float,
        required=True,
        metavar=("MIN", "MAX"),
        help="explicit increasing range for the configured x state variable",
    )
    parser.add_argument(
        "--y-range",
        nargs=2,
        type=float,
        required=True,
        metavar=("MIN", "MAX"),
        help="explicit increasing range for the configured y state variable",
    )
    parser.add_argument(
        "--mesh-points",
        type=int,
        default=201,
        metavar="N",
        help="number of linear sampling points on each axis (default: 201)",
    )
    parser.add_argument(
        "--condition",
        type=_fixed_condition,
        action="append",
        default=[],
        metavar="STATE_VARIABLE=VALUE",
        help="fixed state variable and value; repeat for multiple conditions",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "save the figure to this path; the filename extension selects "
            "the output format"
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="display the figure interactively",
    )
    coloring = parser.add_mutually_exclusive_group()
    coloring.add_argument(
        "--color-component",
        metavar="COMPONENT",
        help=(
            "color by the atomic fraction of this independent component; "
            "defaults to the first independent component"
        ),
    )
    coloring.add_argument(
        "--phase-identity-colors",
        action="store_true",
        help="use discrete stable-phase colors instead of composition",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the generalized phase-diagram plotting command."""
    parser = _parser()
    args = parser.parse_args(argv)
    if args.output is None and not args.show:
        parser.error("at least one of --output or --show is required")

    figure = None
    try:
        fixed_conditions = {}
        for variable, value in args.condition:
            if variable in fixed_conditions:
                parser.error(f"duplicate fixed condition {variable!r}")
            fixed_conditions[variable] = value
        configuration = PhaseDiagramConfiguration.read_json(
            args.configuration
        )
        _, figure, _, _ = _prepare_phase_diagram(
            configuration,
            x_range=tuple(args.x_range),
            y_range=tuple(args.y_range),
            mesh_points=args.mesh_points,
            fixed_conditions=fixed_conditions,
            color_component=args.color_component,
            phase_identity_colors=args.phase_identity_colors,
        )
        if args.output is not None:
            figure.savefig(args.output, bbox_inches="tight")
        if args.show:
            plt.show()
    except (OSError, TypeError, ValueError) as error:
        parser.error(str(error))
    finally:
        if figure is not None:
            plt.close(figure)


if __name__ == "__main__":
    main()
