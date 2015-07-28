
.. Documentation for the header file dolfin/la/GenericPreconditioner.h

.. _programmers_reference_cpp_la_genericpreconditioner:

GenericPreconditioner.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericPreconditioner

    This class provides a common base preconditioners.


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& nullspace)
    
        Set the (approximate) null space of the preconditioner operator
        (matrix). This is required for certain preconditioner types,
        e.g. smoothed aggregation multigrid


    .. cpp:function:: void set_coordinates(const std::vector<double>& x, std::size_t dim)
    
        Set the coordinates of the operator (matrix) rows and geometric
        dimension d. This is can be used by required for certain
        preconditioners, e.g. ML. The input for this function can be
        generated using GenericDofMap::tabulate_all_dofs.


