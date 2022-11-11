class NoInversionSymmetryError(Exception):
    """
    An error class defined for not finding the inversion symmetry center.

    """

    def __str__(self):
        return "The target slab does not have inversion symmetry. \n" \
               "Try to increase the tolerance for symmetry detection."


class SlabOrientationError(Exception):
    """
    An error class defined for wrong slab orientation error.
    """

    def __str__(self):
        return "The slab model provided does not have vacuum in the " \
               "c-direction."


class NonCentralInversionSymmetryError(Exception):
    """
    An error class defined for the slab model which does not have the
    inversion symmetry center located at the origin.
    """

    def __str__(self):
        return "The inversion symmetry center is not at the origin. Please " \
               "check the structure and \n try to move it to the origin."


class PrimitiveStructureFinderError(Exception):
    """
    An error class defined for wrong refined structure (some sites are lost or
    redundant).
    """

    def __str__(self):
        return "The number of sites in the refined structure is not a " \
               "multiple of the number of sites in the structure that" \
               " should be. \n. Please check the primitive structure finder."


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
        """

        """
        return "The input structure is not a slab model. Please double check!"


class PolarSurfaceError(Exception):
    """
    An error class defined for the polar surface.

    Args:
        num_layers: Number of layers in the slab model.
        num_fixed: Number of fixed layers in the slab model.

    """

    def __init__(self, num_layers, num_fixed):
        self.num_layers = num_layers
        self.num_fixed = num_fixed

    def __str__(self):
        if self.num_fixed % 2 == 0:
            return "The polar surface has to have the odd number of fixed" \
                   " layers."
        if (self.num_layers % 2) == 0:
            return "The polar surface has to have the odd number of layers."


class NonPolarSurfaceError(Exception):
    """
    An error class defined for the non-polar surface.

    Args:
        num_layers: Number of layers in the slab model.
        num_fixed: Number of fixed layers in the slab model.

    """

    def __init__(self, num_layers, num_fixed):

        self.num_layers = num_layers
        self.num_fixed = num_fixed

    def __str__(self):
        if (self.num_layers % 2) != 0 and (self.num_fixed % 2 == 0):
            return "The number of fixed layers has to be odd for the " \
                   "non-polar surface with odd number of layers."
        if (self.num_layers % 2) == 0 and (self.num_fixed % 2 != 0):
            return "The number of fixed layers has to be even for the " \
                   "non-polar surface with even number of layers."


class TooLargeSlabError(Exception):
    """
    An error class defined for the too large slab model.
    """

    def __str__(self):
        return "The target cell size defined here is too large which will " \
               "generate too many enumerated structures. \n" \
               "Please consider to decrease the cell size."


class InvalidCompositionError(Exception):
    """
    An error class defined for the invalid user defined composition list to
    be enumerated.
    """

    def __str__(self):
        return "Please double check the composition list or target cell " \
               "size! By applying one of them, \n the num of target " \
               "atom on the surface is not an integer."


class InvalidPhasesAlignError(Exception):
    """
    An error class defined for the invalid phases used to calculate the
    alignment energy.
    """

    def __str__(self):
        return "Please double check the reference phases used here. The " \
               "alignment energy should be \n calculated on the basis of " \
               "two same slab models but with different number of layers."


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
        return "Please double check the symmetric parameter defined in the " \
               "input json file. \n It is not compatible with slab model."


class InvalidInputFormatError(Exception):
    """
    An error class defined to show that the user defined input format is not
    valid.
    """

    def __str__(self):
        return "Please check the input format. \n" \
               "The input should be formatted as JSON file. For details, " \
               "see examples in the example folder."
