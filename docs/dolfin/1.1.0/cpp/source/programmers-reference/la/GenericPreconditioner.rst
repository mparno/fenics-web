
.. Documentation for the header file dolfin/la/GenericPreconditioner.h

.. _programmers_reference_cpp_la_genericpreconditioner:

GenericPreconditioner.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericPreconditioner

    This class provides a common base preconditioners.


    .. cpp:function:: void set_nullspace(const std::vector<const GenericVector*> nullspace)
    
        Set the (approximate) null space of the preconditioner operator
        (matrix). This is required for certain preconditioner types,
        e.g. smoothed aggregation multigrid


