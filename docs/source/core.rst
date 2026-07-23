====
core
====

``SurfaceEnumerator`` class
***************************

.. autoclass:: surface_pd.core.SurfaceEnumerator
    :members:

``SurfaceEnumerationMetadata`` class
************************************

.. autoclass:: surface_pd.core.SurfaceEnumerationMetadata
    :members:

``EnumerationSlab`` class
*************************

The supported slab API consists of the two high-level factories, immutable
analysis, inversion-symmetry inspection, surface-enumeration configuration
properties, and finalized-result metadata. Geometry repair and legacy command
pipeline operations are internal implementation details.

.. autoclass:: surface_pd.core.EnumerationSlab
    :members: from_structure, from_file, analyze, has_inversion_symmetry,
              direction, layer_tolerance_angstrom, enumerated_species,
              num_enumerated_layers, symmetric, enumeration_metadata

``SlabLayer`` class
*******************

.. autoclass:: surface_pd.core.SlabLayer
    :members:

``SlabAnalysis`` class
**********************

.. autoclass:: surface_pd.core.SlabAnalysis
    :members:
