"""
SurfacePD error sub-package

This sub-package consists of a following module:

error.py --> Collection of error class defined by developer to guide user.
"""

from .error import *

__all__ = ['NoInversionSymmetryError', 'SlabOrientationError',
           'NonCentralInversionSymmetryError',
           'PrimitiveStructureFinderError',
           'NonDefinedSelectiveDynamicsError',
           'PolarSurfaceError', 'NonPolarSurfaceError', 'NonSlabError',
           'TooLargeSlabError', 'InvalidCompositionError']
