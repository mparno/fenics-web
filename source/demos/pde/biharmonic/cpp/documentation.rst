.. Documentation for the biharmonic demo from DOLFIN.

.. _demos_cpp_pde_biharmonic:


Biharmonic equation
===================

.. include:: ../common.txt


Implementation
--------------

The implementation is split in two files, a form file containing the definition
of the variational forms expressed in UFL and the solver which is implemented
in a C++ file.

Complete code
-------------

Complete UFL file
^^^^^^^^^^^^^^^^^

.. literalinclude:: Biharmonic.ufl
   :start-after: # Compile
   :language: python

Complete main file
^^^^^^^^^^^^^^^^^^

.. literalinclude:: main.cpp
   :start-after: // using
   :language: c++
