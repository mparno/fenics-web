.. Documentation for the header file dolfin/la/DefaultFactory.h

.. _programmers_reference_cpp_la_defaultfactory:

DefaultFactory.h
================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: DefaultFactory

    *Parent class*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: DefaultFactory()
    
        Constructor

    .. cpp:function:: virtual dolfin::GenericLinearSolver*
                                                           create_krylov_solver(std::string method, std::string pc) const
    
        Create Krylov solver

    .. cpp:function:: virtual dolfin::GenericLinearSolver* create_lu_solver() const
    
        Create LU solver

    .. cpp:function:: virtual dolfin::GenericMatrix* create_matrix() const
    
        Create empty matrix

    .. cpp:function:: virtual dolfin::GenericSparsityPattern* create_pattern() const
    
        Create empty sparsity pattern

    .. cpp:function:: virtual dolfin::GenericVector* create_local_vector() const
    
        Create empty vector (local)

    .. cpp:function:: virtual dolfin::GenericVector* create_vector() const
    
        Create empty vector (global)

    .. cpp:function:: virtual ~DefaultFactory()
    
        Destructor

