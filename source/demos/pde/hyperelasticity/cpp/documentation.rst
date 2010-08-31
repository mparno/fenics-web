.. Documentation for the hyperelasticity demo from DOLFIN.

.. _demos_cpp_pde_hyperelasticity:

Hyperelasticity
===============

.. include:: ../common.txt

Implementation
--------------

The implementation is split in two files, a form file containing the
definition of the variational forms expressed in UFL and the solver
which is implemented in a C++ file.


Complete code
-------------

.. literalinclude:: main.cpp
   :start-after: // Begin demo
   :language: c++

