.. _installation:
.. index:: Installation


Installation
**************

Download
--------------

To download :program:`PyCatKin`, you can clone the Git repo::

	git clone https://github.com/aab64/PyCatKin.git

or download a .zip archive::

    https://github.com/aab64/PyCatKin/archive/refs/heads/master.zip

Setup
------

:program:`PyCatKin` can be installed by running::

	pip install .

from the root directory.

Dependencies
--------------

:program:`PyCatKin` depends on the following `Python` packages:

 - `NumPy <https://www.numpy.org/>`_ used for array manipulations.
 - `SciPy <https://www.scipy.org/>`_ used for integration and root finding.
 - `matplotlib <https://matplotlib.org/>`_ for plotting results.
 - `ASE <https://wiki.fysik.dtu.dk/ase>`_ for handling atomic structures.
 - `pandas <https://pandas.pydata.org/>`_ for handling file output.
 - `graphviz <https://graphviz.org/>`_ (optional) for profiling the code, must be added to the path.
