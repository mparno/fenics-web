.. Documentation for the header file dolfin/la/LinearAlgebraFactory.h

.. _programmers_reference_cpp_la_linearalgebrafactory:

LinearAlgebraFactory.h
======================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: LinearAlgebraFactory

    .. cpp:function:: LinearAlgebraFactory()
    
        Constructor

    .. cpp:function:: dolfin::GenericMatrix* create_matrix() const = 0
    
        Create empty matrix

    .. cpp:function:: dolfin::GenericVector* create_vector() const = 0
    
        Create empty vector (global)

    .. cpp:function:: dolfin::GenericVector* create_local_vector() const = 0
    
        Create empty vector (local)

    .. cpp:function:: dolfin::GenericSparsityPattern* create_pattern() const = 0
    
        Create empty sparsity pattern (returning zero if not used/needed)

    .. cpp:function:: dolfin::GenericLinearSolver* create_lu_solver() const = 0
    
        Create LU solver

    .. cpp:function:: dolfin::GenericLinearSolver* create_krylov_solver(std::string method, std::string pc) const = 0
    
        Create Krylov solver

