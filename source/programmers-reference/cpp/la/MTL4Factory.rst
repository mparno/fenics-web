.. Documentation for the header file dolfin/la/MTL4Factory.h

.. _programmers_reference_cpp_la_mtl4factory:

MTL4Factory.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MTL4Factory

    *Parent class*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: GenericSparsityPattern* create_pattern() const
    
        Dummy sparsity pattern

    .. cpp:function:: ITLKrylovSolver* create_krylov_solver(std::string method,
                                                            std::string pc) const
    
        Create Krylov solver

    .. cpp:function:: MTL4Matrix* create_matrix() const
    
        Create empty matrix

    .. cpp:function:: MTL4Vector* create_local_vector() const
    
        Create empty vector (local)

    .. cpp:function:: MTL4Vector* create_vector() const
    
        Create empty vector (global)

    .. cpp:function:: UmfpackLUSolver* create_lu_solver() const
    
        Create LU solver

