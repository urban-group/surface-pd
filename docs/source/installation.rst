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
    the ``enum.x`` and ``makestr.x`` must be in the path. The detailed
    compilation steps are fully described at https://github.com/msg-byu/enumlib.

    The compilation step can also make use of the pymatgen command line tool as
    follows (assume that the pymatgen has been installed successfully): ::

        $ pmg config --install enumlib

Then put these in your PATH somewhere.

Installation test (This needs to further check)
***********************************************

To check whether the **surface-pd** package is successfully installed, you
can import **surface-pd** within Python interpreter::

    >>> import surface-pd

Or the following command can be executed (with an error in the unprepared state)::

    $ surface-pd


