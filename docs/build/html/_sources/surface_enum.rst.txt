.. _label_surface_enumeration:

===================
Surface enumeration
===================


Workflow
********

.. image:: images/flowchart.png
    :width: 500
    :align: center

- The compound, shape, geometry, sites property of the input slab model will be checked first to make sure the user defined slab model is valid.
- The to-be-enumerated composition that user defined will be further checked to make sure it will yield integral number of atoms.
- The slab model will be separated as the top, center, and bottom regions according to the statement of constrains defined for all atoms (user defined).
- The top and center regions will be combined and create a pseudo slab model.
- The target atoms in the pseudo slab model will be substituted by dummy species to facilitate the EnumWithComposition class to detect and enumerate.
- The "top" surface of the pseudo slab model will be enumerated in 3D to generate all geometrically distinct slab models (a, b, and c lattice parameters are all possible to change). 
- Applying the filtration steps to remove those enumerated pseudo slab models have c lattice parameter changed since we only want 2D surface enumeration (only a and b lattice parameters are allowed to change).
- The remaining enumerated pseudo slab models will be symmetrized on the basis of the inversion symmetry.
- The symmetrized "real" slab model will further be refined, which means the atoms will be moved to the expected symmetry positions.
- The inversion symmetry center will be shifted to the origin (0, 0, 0) and let VASP determines the point group symmetry and the space group.

Enumlib code
************

The enumlib code was developed by Prof. Gus Hart's group by Brigham Young
University (BYU).

For the full description of the enumlib code, please refer to its
`github link <https://github.com/msg-byu/enumlib>`__ and the following papers.

1. `Algorithm for generating derivative structures
<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.77.224115>`__

2. `Generating derivative structures from multilattices: Algorithm and application to hcp alloys
<https://journals.aps.org/prb/abstract/10.1103/PhysRevB.80.014120>`__

3. `Generating derivative structures at a fixed concentration
<https://www.sciencedirect.com/science/article/abs/pii/S092702561200081X>`__

4. `Generating derivative superstructures for systems with high configurational freedom
<https://www.sciencedirect.com/science/article/abs/pii/S0927025617302069>`__
