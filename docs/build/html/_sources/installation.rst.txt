============
Installation
============

Requirements
************
* `Monty <https://pythonhosted.org/monty/index.html>`_
* `Python>=3 <https://www.python.org/>`_
* `Numpy <https://numpy.org/doc/stable/>`_
* `Matplotlib <https://matplotlib.org/>`_
* `Pandas <https://pandas.pydata.org/>`_
* `Pymatgen <https://pymatgen.org/>`_
* `Enumlib <https://github.com/msg-byu/enumlib/>`_

Installation using pip
**********************

**surface-pd** is a Python 3 package and the stable version can be installed
via pip_. ::

    $ pip install ***

.. _PIP: https://pip.pypa.io/en/stable/

Installation from source code
*****************************

.. :Git clone:

The package can also be built from `surface-pd
GitHub <https://github.com/urban-group/surface-pd>`_ by git
clone
and setup::

    $ git clone git@github.com:urban-group/surface-pd.git
    $ cd surface-pd
    $ pip install --user .

Additional requirements
***********************

.. note::

    In order to use the enumeration functionalities provided in this package,
    the ``enum.x`` and ``makestr.x`` must be in the path.

    You can run the following several lines of code to compile the enumlib
    source code in your terminal. ::

        $ git clone  --branch v2.0.4 --recursive https://github.com/msg-byu/enumlib.git
        $ cd enumlib/symlib/src && make F90=gfortran
        $ cd enumlib/src && make F90=gfortran
        $ cd enumlib/src && F90=gfortran make enum.x
        $ cd enumlib/src && F90=gfortran make makestr.x


    After the compilation, you should add these two executable files (``enum.x`` and ``makestr.x``) to your
    global path . ::

        $ cp enumlib/src/*.x "your_global_path"

    To check whether the enumeration utilities are available: ::

        $ which enum.x makestr.x

    The paths to these two files should be located.



Installation test
***********************************************

After installation, you should be able to see three scripts, they
are::

    $ surface-enumeration.py
    $ surface-pd-plot.py
    $ discharge_pd_gene.py

To check whether the **surface-pd** package is successfully installed, you
can try if the following command can be executed::

    $ surface-enumeration.py -h

The user-friendly command-line interfaces with help and usage messages
should be automatically generated. Something like::

    $  surface-enumeration.py -h
    usage: surface-enumeration.py [-h] [--max-structures MAX_STRUCTURES] [--generate-poscar] json_file_path
    This code will enumerate the input slab model with user defined target
    species and composition. The input slab structure should be as small as
    possible because the unit cell of the slab will also be enumerated, and it
    will increase the number of possibilities that the enumeration can have.
    For the detailed algorithm behind the enumeration, please see the
    following references below.
    (1) Morgan, W. S.; Hart, G. L. W.; Forcade, R. W.
    Computational Materials Science 2017, 136, 144–149.
    https://doi.org/10.1016/j.commatsci.2017.04.015.
    ...
    ...

