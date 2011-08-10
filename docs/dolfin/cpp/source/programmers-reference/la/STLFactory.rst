
.. Documentation for the header file dolfin/la/STLFactory.h

.. _programmers_reference_cpp_la_stlfactory:

STLFactory.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: STLFactory

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: STLMatrix* create_matrix() const
    
        Create empty matrix


    .. cpp:function:: uBLASVector* create_vector() const
    
        Create empty vector (global)


    .. cpp:function:: uBLASVector* create_local_vector() const
    
        Create empty vector (local)


    .. cpp:function:: GenericSparsityPattern* create_pattern() const
    
        Create empty sparsity pattern


    .. cpp:function:: GenericLinearSolver* create_lu_solver() const
    
        Create LU solver


    .. cpp:function:: GenericLinearSolver* create_krylov_solver(std::string method, std::string pc) const
    
        Create Krylov solver


    .. cpp:function:: static STLFactory& instance()
    
        Return singleton instance


    .. cpp:function:: STLFactory()
    
        Private Constructor


