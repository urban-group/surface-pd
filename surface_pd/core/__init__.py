"""
Core modules for surface-pd slab enumeration.

This subpackage contains the core functionality for creating, checking, and
enumerating slab models for surface phase diagram generation.

Modules
-------
enumeration_slab
    EnumerationSlab class for surface enumeration structures.
enum
    Enumeration classes for systematic surface composition generation.
pre_check
    Pre-processing validation for input slab structures.
post_check
    Post-processing validation for enumerated structures.
"""

from .enum import EnumWithComposition
from .enumeration_slab import EnumerationSlab, SlabAnalysis, SlabLayer

__all__ = [
    "EnumerationSlab",
    "SlabAnalysis",
    "SlabLayer",
    "EnumWithComposition",
]
