.. Documentation for the header file dolfin/la/LinearAlgebraFactory.h

.. _programmers_reference_cpp_la_Mesh:

LinearAlgebraFactory.h
======================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: LinearAlgebraFactory

    .. cpp:function:: LinearAlgebraFactory()
    
        Constructor

    .. cpp:function:: virtual dolfin::GenericLinearSolver*
                                                           create_krylov_solver(std::string method, std::string pc) const = 0
    
        Create Krylov solver

    .. cpp:function:: virtual dolfin::GenericLinearSolver* create_lu_solver() const = 0
    
        Create LU solver

    .. cpp:function:: virtual dolfin::GenericMatrix* create_matrix() const = 0
    
        Create empty matrix

    .. cpp:function:: virtual dolfin::GenericSparsityPattern* create_pattern() const = 0
    
        Create empty sparsity pattern (returning zero if not used/needed)

    .. cpp:function:: virtual dolfin::GenericVector* create_local_vector() const = 0
    
        Create empty vector (local)

    .. cpp:function:: virtual dolfin::GenericVector* create_vector() const = 0
    
        Create empty vector (global)

    .. cpp:function:: virtual ~LinearAlgebraFactory()
    
        Destructor

