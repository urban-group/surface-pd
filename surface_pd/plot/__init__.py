"""
Plotting subpackage for surface phase diagram generation and visualization.

This subpackage renders evaluated phase-diagram results using explicit,
chemistry-independent axes and optional composition coloring.
"""

from .phase_diagram import CompositionColoring, plot_phase_diagram

__all__ = [
    "CompositionColoring",
    "plot_phase_diagram",
]
