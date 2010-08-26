.. Documentation for the header file dolfin/la/uBLASFactory.h

.. _programmers_reference_cpp_la_Mesh:

uBLASFactory.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: uBLASFactory

    *Parent class*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: GenericLinearSolver* create_krylov_solver(std::string method, std::string pc) const
                                                                //
    
        Create Krylov solver

    .. cpp:function:: SparsityPattern* create_pattern() const
    
        Create empty sparsity pattern

    .. cpp:function:: UmfpackLUSolver* create_lu_solver() const
    
        Create LU solver

    .. cpp:function:: uBLASMatrix<Mat>* create_matrix() const
    
        Create empty matrix

    .. cpp:function:: uBLASVector* create_local_vector() const
    
        Create empty vector (local)

    .. cpp:function:: uBLASVector* create_vector() const
    
        Create empty vector

    .. cpp:function:: virtual ~uBLASFactory()
    
        Destructor

