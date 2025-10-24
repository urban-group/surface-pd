"""
Custom error classes for surface-pd package.

This module defines specific error types for different failure modes in
surface enumeration and analysis, providing clear error messages to guide
users in resolving issues.
"""


class NoInversionSymmetryError(Exception):
    """An error class defined for not finding the inversion symmetry center."""

    def __str__(self):
        return (
            "The target slab does not have inversion symmetry. \n"
            "Try to increase the tolerance for symmetry detection."
        )


class SlabOrientationError(Exception):
    """An error class defined for wrong slab orientation error."""

    def __str__(self):
        return (
            "The slab model given by user does not have vacuum in the "
            "c-direction."
        )


class NonCentralInversionSymmetryError(Exception):
    """
    An error class defined for the slab model which does not have the
    inversion symmetry center located at the origin.
    """

    def __str__(self):
        return (
            "The inversion symmetry center is not at the origin. Please "
            "check the structure and \n try to move it to the origin."
        )


class PrimitiveStructureFinderError(Exception):
    """
    An error class defined for wrong refined structure (some sites are lost or
    redundant).
    """

    def __str__(self):
        return (
            "The number of sites in the refined structure is not a "
            "multiple of the number of sites in the structure that"
            "is going to be enumerated."
        )


class NonDefinedSelectiveDynamicsError(Exception):
    """
    An error class defined for the slab model that does not have selective
    dynamics defined for all sites (selective dynamics is extremely useful
    in several steps).
    """

    def __str__(self):
        return "Not all sites in the slab model have selective dynamics."


class NonSlabError(Exception):
    """
    An error class defined for the structure that is not a slab model
    (no enough vacuum region in the model).
    """

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
    """An error class defined for the too large slab model."""

    def __str__(self):
        return (
            "The target cell size defined here is too large which will "
            "generate too many enumerated structures. \n"
            "Please decrease the target cell size."
        )


class InvalidCompositionError(Exception):
    """
    An error class defined for the invalid user defined composition list to
    be enumerated.
    """

    def __str__(self):
        return (
            "Please double check the composition list or target cell "
            "size! By applying one of them, \n the num of target "
            "species should be an integer."
        )


class InvalidPhasesAlignError(Exception):
    """
    An error class defined for the invalid phases used to calculate the
    alignment energy.
    """

    def __str__(self):
        return (
            "Please double check the reference phases used here. The "
            "alignment energy should be \n calculated on the basis of "
            "same amount of surface area (same supercell size)."
        )


class NonIntegerError(Exception):
    """
    An error class defined to check whether the after enumerated number of
    atoms is an integer with tolerance.
    """

    def __str__(self):
        return "The after enumerated number of atoms is not an integer."


class IncompatibleSymmError(Exception):
    """
    An error class defined to show that the user defined symmetric parameter is
     not compatible with code determined symmetric parameter.
    """

    def __str__(self):
        return (
            "Please double check the symmetric parameter defined in the "
            "input json file. \n It is not compatible with slab model."
        )


class InvalidInputFormatError(Exception):
    """
    An error class defined to show that the user defined input format is not
    valid.
    """

    def __str__(self):
        return (
            "Please check the input format. \n"
            "The input should be formatted as JSON file. For details, "
            "please check the documentation or the examples."
        )
