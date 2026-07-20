"""
SurfacePD error sub-package.

This sub-package consists of a following module:

error.py --> Collection of error class defined by developer to guide user.
"""

from .error import (
    IncompatibleSymmError,
    InvalidCompositionError,
    InvalidInputFormatError,
    NoInversionSymmetryError,
    NonCentralInversionSymmetryError,
    NonDefinedSelectiveDynamicsError,
    NonIntegerError,
    NonSlabError,
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
    "NonSlabError",
    "TooLargeSlabError",
    "InvalidCompositionError",
    "IncompatibleSymmError",
    "NonIntegerError",
    "InvalidInputFormatError",
]
