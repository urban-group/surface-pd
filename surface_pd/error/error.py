class NoInversionSymmetryError(Exception):
    """

    """
    def __str__(self):
        """
        Error message for not finding the inversion symmetry center
        Returns:

        """
        return "The target slab does not have inversion symmetry. \n" \
               "Try to increase the tolerance for symmetry detection."

class SlabOrientationError(Exception):
    """

    """
    def __str__(self):
        return "The slab model provided is not a tall cuboid. Plaese " \
               "consider to construct a slab with c > (a and b)."

class NonCentralInversionSymmetryError(Exception):
    """

    """
    def __str__(self):
        return "The inversion symmetry center is not at the origin. Please " \
               "check the structure and \n try to move it to the origin."


class PrimitiveStructureFinderError(Exception):
    """

    """
    def __str__(self):
        return "The number of sites in the refined structure is not a " \
               "multiple of the number of sites in the structure that" \
               " should be. \n. Please check the primitive structure finder."


class NonDefinedSelectiveDynamicsError(Exception):
    """

    """
    def __str__(self):
        return "Not all sites in the slab model have selective dynamics."


class NonSlabError(Exception):
    """

    """
    def __str__(self):
        """

        """
        return "The input structure is not a slab model. Please double check!"


class PolarSurfaceError(Exception):
    """

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

    """
    def __str__(self):
        return "The target cell size defined here is too large which will " \
               "generate too many enumerated structures. \n" \
               "Please consider to decrease the cell size."