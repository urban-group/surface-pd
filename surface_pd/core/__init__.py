"""
Core modules for surface-pd slab enumeration.

This subpackage contains the core functionality for creating, checking, and
enumerating slab models for surface phase diagram generation.

Modules
-------
slab
    Slab class for representing and manipulating surface slab structures.
enum
    Enumeration classes for systematic surface composition generation.
pre_check
    Pre-processing validation for input slab structures.
post_check
    Post-processing validation for enumerated structures.
"""

from .enum import EnumWithComposition
from .post_check import PostCheck
from .pre_check import PreCheck
from .slab import Slab

__all__ = ["Slab", "EnumWithComposition", "PreCheck", "PostCheck"]
