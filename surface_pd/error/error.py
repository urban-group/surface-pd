"""
Custom error classes for surface-pd package.

This module defines specific error types for different failure modes in
surface enumeration and analysis, providing clear error messages to guide
users in resolving issues.
"""


class NoInversionSymmetryError(Exception):
    """Raised when an inversion symmetry center cannot be found."""

    def __str__(self):
        return (
            "The target slab does not have inversion symmetry. \n"
            "Try to increase the tolerance for symmetry detection."
        )


class SlabOrientationError(Exception):
    """Raised when slab vacuum is not aligned with the expected direction."""

    def __str__(self):
        return (
            "The slab model given by user does not have vacuum in the "
            "c-direction."
        )


class NonCentralInversionSymmetryError(Exception):
    """Raised when a slab's inversion center is not at the origin."""

    def __str__(self):
        return (
            "The inversion symmetry center is not at the origin. Please "
            "check the structure and \n try to move it to the origin."
        )


class PrimitiveStructureFinderError(Exception):
    """Raised when refinement loses sites or introduces redundant sites."""

    def __str__(self):
        return (
            "The number of sites in the refined structure is not a "
            "multiple of the number of sites in the structure that"
            "is going to be enumerated."
        )


class NonDefinedSelectiveDynamicsError(Exception):
    """Raised when some slab sites lack selective-dynamics flags."""

    def __str__(self):
        return "Not all sites in the slab model have selective dynamics."


class NonSlabError(Exception):
    """Raised when a structure lacks enough vacuum to be a slab model."""

    def __str__(self):
        """Return error message."""
        return "The input structure is not a slab model. Please double check!"


class PolarSurfaceError(Exception):
    """An error class defined for the polar surface error."""

    def __init__(self, num_layers, num_fixed):
        self.num_layers = num_layers
        self.num_fixed = num_fixed

    def __str__(self):
        if self.num_fixed % 2 == 0:
            return (
                "The polar surface has to have the odd number of fixed layers."
            )
        else:
            return (
                f"There are {self.num_layers} layers, "
                f"with {self.num_fixed} fixed layers."
            )


class NonPolarSurfaceError(Exception):
    """An error class defined for the non-polar surface error."""

    def __init__(self, num_layers, num_fixed):
        self.num_layers = num_layers
        self.num_fixed = num_fixed

    def __str__(self):
        if (self.num_layers % 2) != 0 and (self.num_fixed % 2 == 0):
            return (
                "The number of fixed layers has to be odd for the "
                "non-polar surface."
            )
        else:
            return (
                f"There are {self.num_layers} layers, "
                f"with {self.num_fixed} fixed layers."
            )


class TooLargeSlabError(Exception):
    """Raised when a requested cell would produce too many structures."""

    def __str__(self):
        return (
            "The target cell size defined here is too large which will "
            "generate too many enumerated structures. \n"
            "Please decrease the target cell size."
        )


class InvalidCompositionError(Exception):
    """Raised when a composition is incompatible with the target cell size."""

    def __str__(self):
        return (
            "Please double check the composition list or target cell "
            "size! By applying one of them, \n the num of target "
            "species should be an integer."
        )


class NonIntegerError(Exception):
    """Raised when an enumerated atom count is not effectively integral."""

    def __str__(self):
        return "The after enumerated number of atoms is not an integer."


class IncompatibleSymmError(Exception):
    """Raised when requested symmetry conflicts with the slab geometry."""

    def __str__(self):
        return (
            "Please double check the symmetric parameter defined in the "
            "input json file. \n It is not compatible with slab model."
        )


class InvalidInputFormatError(Exception):
    """Raised when CLI input cannot be parsed as the expected JSON format."""

    def __str__(self):
        return (
            "Please check the input format. \n"
            "The input should be formatted as JSON file. For details, "
            "please check the documentation or the examples."
        )
