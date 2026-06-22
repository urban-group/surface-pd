"""
SurfacePD error sub-package.

This sub-package consists of a following module:

error.py --> Collection of error class defined by developer to guide user.
"""

from .error import (
    IncompatibleSymmError,
    InvalidCompositionError,
    InvalidInputFormatError,
    InvalidPhasesAlignError,
    NoInversionSymmetryError,
    NonCentralInversionSymmetryError,
    NonDefinedSelectiveDynamicsError,
    NonIntegerError,
    NonPolarSurfaceError,
    NonSlabError,
    PolarSurfaceError,
    PrimitiveStructureFinderError,
    SlabOrientationError,
    TooLargeSlabError,
)

__all__ = [
    "NoInversionSymmetryError",
    "SlabOrientationError",
    "NonCentralInversionSymmetryError",
    "PrimitiveStructureFinderError",
    "NonDefinedSelectiveDynamicsError",
    "PolarSurfaceError",
    "NonPolarSurfaceError",
    "NonSlabError",
    "TooLargeSlabError",
    "InvalidCompositionError",
    "InvalidPhasesAlignError",
    "IncompatibleSymmError",
    "NonIntegerError",
    "InvalidInputFormatError",
]
