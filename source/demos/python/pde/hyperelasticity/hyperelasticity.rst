.. Documentation for the hyperelasticity demo from DOLFIN.

.. _demos_python_pde_hyperelasticity:

Hyperelasticity
===============

.. include:: ../../../common/pde/hyperelasticity/hyperelasticity.txt


Implementation
--------------

This demo is implemented in a single Python file, :download:`demo.py`, which
contains both the variational forms and the solver.

First, the ``dolfin`` module is imported:

.. code-block:: python

    from dolfin import *
