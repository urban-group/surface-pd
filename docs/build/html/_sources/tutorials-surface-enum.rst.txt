==============================
Slab Model Surface Enumeration
==============================


Slab model preparation
**********************
The slab model can be truncated from the bulk structure model. Assume now,
you already have a slab model. The DFT structure input file in VASP format
should like this. ::

    Li4 Co3 O8
    1.0
    2.832550 0.000000 0.000000
    -1.416275 2.453061 0.000000
    0.000000 0.000000 31.000000
    Li Co O
    4 3 8
    Selective dynamics
    direct
    0.000000 0.000000 0.199664 F F F Li # atom index = 0
    0.333333 0.666667 0.047755 T T T Li # atom index = 1
    0.666667 0.333333 0.351573 F F F Li # atom index = 2
    0.333333 0.666667 0.503483 T T T Li # atom index = 3
    0.666667 0.333333 0.123709 T T T Co # atom index = 4
    0.333333 0.666667 0.275619 F F F Co # atom index = 5
    0.000000 0.000000 0.427529 T T T Co # atom index = 6
    0.333333 0.666667 0.156909 T T T O  # atom index = 7
    0.000000 0.000000 0.090509 T T T O  # atom index = 8
    0.666667 0.333333 0.242419 F F F O  # atom index = 9
    0.666667 0.333333 0.005000 T T T O  # atom index = 10
    0.000000 0.000000 0.546237 T T T O  # atom index = 11
    0.000000 0.000000 0.308819 F F F O  # atom index = 12
    0.666667 0.333333 0.460728 T T T O  # atom index = 13
    0.333333 0.666667 0.394328 T T T O  # atom index = 14

Note that the surface regions are defined by having the selective dynamics
labeled at the end, represented by "T T T" as relaxed and "F F F" as fixed.
These properties are required to build a slab model.

.. note::
    Please make sure that
    the input slab model has the statement of constraints because the code takes
    use of the constraints to distinguish relaxed and fixed regions in the slab
    model.

.. |LiCoO2| replace:: LiCoO\ :sub:`2`

The following figures show the example (001) and (104) surface slab models
of |LiCoO2|.

.. image:: images/slab-model-demo.png
    :width: 400
    :align: center

\

Command line execution
**********************

The general usage of the code can be seen by calling: ::

    $ surface-enumeration.py --help

All of the available options and their default values will be shown.

General long format::

    $ surface-enumeration.py 'slab model' --target-cell-size 'int' --lithium-composition 'float(0<=x<=1)' --oxygen-composition 'float(0<=y<=1)' --generate-poscar



General short format::

    $ surface-enumeration.py 'slab model' -s 'int' -L 'float(0<=x<=1)' -O 'float(0<=y<=1)' -g

For example, you can go to ``surface-pd/example/enumeration-examples/LCO``
and try to do a full enumeration on the (104) surface slab model of
|LiCoO2|. ::

    $ surface-enumeration.py LCO_terminated_LCO_10L15A.vasp -s 2

For this slab model, it has 2 Li and O atoms in the surface
relaxed region. Since the target number of Li and O atoms is set to 4 for
this example case, you can always see the ``--target-cell-size`` option
turns on and is set to 2 in the following examples in order to create the
supercell.

If everything goes well, you should see the following: ::

    target_cell_size = 2
    Composition of lithium on the surface will be [1.0, 0.75, 0.5, 0.25, 0.0].
    Composition of oxygen on the surface will be [1.0, 0.75, 0.5, 0.25, 0.0].
    Scaling matrix used here is: [2, 1, 1]
    The enumeration found 1 (0+1) distinct structures for 100.0% Li and 100.0% O.
    The enumeration found 4 (4+0) distinct structures for 100.0% Li and 75.0% O.
    The enumeration found 6 (6+0) distinct structures for 100.0% Li and 50.0% O.
    The enumeration found 4 (4+0) distinct structures for 100.0% Li and 25.0% O.
    The enumeration found 1 (0+1) distinct structures for 100.0% Li and 0.0% O.
    The enumeration found 4 (4+0) distinct structures for 75.0% Li and 100.0% O.
    The enumeration found 16 (16+0) distinct structures for 75.0% Li and 75.0% O.
    The enumeration found 24 (24+0) distinct structures for 75.0% Li and 50.0% O.
    The enumeration found 16 (16+0) distinct structures for 75.0% Li and 25.0% O.
    The enumeration found 4 (0+4) distinct structures for 75.0% Li and 0.0% O.
    The enumeration found 6 (6+0) distinct structures for 50.0% Li and 100.0% O.
    The enumeration found 24 (24+0) distinct structures for 50.0% Li and 75.0% O.
    The enumeration found 36 (36+0) distinct structures for 50.0% Li and 50.0% O.
    The enumeration found 24 (24+0) distinct structures for 50.0% Li and 25.0% O.
    The enumeration found 6 (0+6) distinct structures for 50.0% Li and 0.0% O.
    The enumeration found 4 (4+0) distinct structures for 25.0% Li and 100.0% O.
    The enumeration found 16 (16+0) distinct structures for 25.0% Li and 75.0% O.
    The enumeration found 24 (24+0) distinct structures for 25.0% Li and 50.0% O.
    The enumeration found 16 (16+0) distinct structures for 25.0% Li and 25.0% O.
    The enumeration found 4 (0+4) distinct structures for 25.0% Li and 0.0% O.
    The enumeration found 1 (0+1) distinct structures for 0.0% Li and 100.0% O.
    The enumeration found 4 (4+0) distinct structures for 0.0% Li and 75.0% O.
    The enumeration found 6 (6+0) distinct structures for 0.0% Li and 50.0% O.
    The enumeration found 4 (4+0) distinct structures for 0.0% Li and 25.0% O.
    The enumeration found 1 (0+1) distinct structures for 0.0% Li and 0.0% O.
    256 distinct structures are found totally.

If you only want to enumerate the slab model with 50% coverage of Li and O
atoms on the surface relaxed region, you can try: ::

    $ surface-enumeration.py LCO_terminated_LCO_10L15A.vasp -s 2 -L 0.5 -O 0.5

If the ``--generate-poscar`` optional argument is defined, you should be
able to see the saved enumerated slab models locally. The slab models are
stored in **VESTA** format and should be able visualize via
`VESTA <https://jp-minerals.org/vesta/en/>`__.

