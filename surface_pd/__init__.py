"""
surface-pd: Surface Phase Diagram Generator.

An open-source Python package designed to automate surface reconstruction
enumeration and plotting of surface phase diagrams for materials science
applications.

This package uses density functional theory (DFT) energies to calculate the
surface free energy of each reconstruction. It can be directly used for
surface enumeration of catalysts, electrodes, and metal surfaces, as well
as for studying corrosion on surfaces.

Modules
-------
analysis
    Analysis modules to analyze and manipulate slab models.
core
    Core modules to check, build, enumerate, and symmetrize slab models.
error
    Custom error classes for quick problem identification.
plot
    Plotting modules to manage DFT results and construct surface phase
    diagrams.
util
    Utility functions for various operations.

The supported domain APIs are available from :mod:`surface_pd.core`,
:mod:`surface_pd.plot`, and :mod:`surface_pd.error`. See the documentation at
https://surface-pd.readthedocs.io for complete examples.

Authors
-------
Xinhao Li (xinhao.li@columbia.edu)
Alexander Urban (a.urban@columbia.edu)
"""

__version__ = "1.0.0"
__all__ = ["__version__"]
