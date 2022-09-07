============
Installation
============

Requirements
************

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
    $ python setup.py install

Additional requirements
***********************

.. note::

    In order to use the enumeration functionalities provided in this package,
    the ``enum.x`` and ``makestr.x`` must be in the path.

    You can run the following several lines of code to complie the enumlib source code in your termnal. ::

        $ git clone  --branch v2.0.4 --recursive https://github.com/msg-byu/enumlib.git
        $ cd enumlib/symlib/src && make F90=gfortran
        $ cd enumlib/src && make F90=gfortran
        $ cd enumlib/src && F90=gfortran make enum.x
        $ cd enumlib/src && F90=gfortran make makestr.x
        $ mkdir -p /opt/bin
        $ cp enumlib/src/*.x /opt/bin 

    The last line of code is just to copy the excecutable ``enum.x`` and ``makestr.x`` files to the path that we think should be global. But you can dynamically change it to anywhere you want.
    The detailed compilation steps are fully described at https://github.com/msg-byu/enumlib.

    To check whether the enumeration utilities are available: ::

        $ which enum.x makestr.x

    The paths to these two files should be located.



Installation test (This needs to further check)
***********************************************

To check whether the **surface-pd** package is successfully installed, you
can import **surface-pd** within Python interpreter::

    >>> import surface-pd

Or the following command can be executed (with an error in the unprepared state)::

    $ surface-pd


