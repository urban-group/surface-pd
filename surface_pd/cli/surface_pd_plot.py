#!/usr/bin/env python3

"""Plot a generalized phase diagram from versioned JSON configuration."""

import argparse
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.figure import Figure

from surface_pd.configuration import PhaseDiagramConfiguration
from surface_pd.plot import plot_phase_diagram
from surface_pd.thermodynamics import PhaseDiagramResult


def _prepare_phase_diagram(
    configuration: PhaseDiagramConfiguration,
) -> tuple[PhaseDiagramResult, Figure, Axes, Colorbar]:
    """Evaluate and render without displaying or saving the figure."""
    if not isinstance(configuration, PhaseDiagramConfiguration):
        raise TypeError(
            "configuration must be a PhaseDiagramConfiguration"
        )
    datasets = configuration.load_datasets()
    result = configuration.diagram_specification.evaluate(
        configuration.model,
        datasets,
    )
    figure, axes, colorbar = plot_phase_diagram(
        result,
        coloring=configuration.create_coloring(),
        cmap=configuration.colormap,
        invert_x_axis=configuration.invert_x_axis,
        invert_y_axis=configuration.invert_y_axis,
    )
    return result, figure, axes, colorbar


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
    return parser


def main(argv: list[str] | None = None) -> None:
    """Run the generalized phase-diagram plotting command."""
    parser = _parser()
    args = parser.parse_args(argv)
    if args.output is None and not args.show:
        parser.error("at least one of --output or --show is required")

    figure = None
    try:
        configuration = PhaseDiagramConfiguration.read_json(
            args.configuration
        )
        _, figure, _, _ = _prepare_phase_diagram(configuration)
        if args.output is not None:
            figure.savefig(args.output)
        if args.show:
            plt.show()
    except (OSError, TypeError, ValueError) as error:
        parser.error(str(error))
    finally:
        if figure is not None:
            plt.close(figure)


if __name__ == "__main__":
    main()
