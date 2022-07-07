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
        return "The slab model provided is not a tall cuboid. Please " \
               "consider to construct a slab with c > (a and b)."


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
    dynamics defined for all sites.
    """

    def __str__(self):
        return "Not all sites in the slab model have selective dynamics."


class NonSlabError(Exception):
    """
    An error class defined for the structure that is not a slab model.
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
