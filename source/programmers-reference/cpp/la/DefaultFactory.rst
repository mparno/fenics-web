.. Documentation for the header file dolfin/la/DefaultFactory.h

.. _programmers_reference_cpp_la_defaultfactory:

DefaultFactory.h
================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: DefaultFactory

    *Parent class*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: DefaultFactory()
    
        Constructor

    .. cpp:function:: dolfin::GenericLinearSolver* create_krylov_solver(std::string method, std::string pc) const
    
        Create Krylov solver

    .. cpp:function:: dolfin::GenericLinearSolver* create_lu_solver() const
    
        Create LU solver

    .. cpp:function:: dolfin::GenericMatrix* create_matrix() const
    
        Create empty matrix

    .. cpp:function:: dolfin::GenericSparsityPattern* create_pattern() const
    
        Create empty sparsity pattern

    .. cpp:function:: dolfin::GenericVector* create_local_vector() const
    
        Create empty vector (local)

    .. cpp:function:: dolfin::GenericVector* create_vector() const
    
        Create empty vector (global)

