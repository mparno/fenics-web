
.. Documentation for the header file dolfin/la/uBLASFactory.h

.. _programmers_reference_cpp_la_ublasfactory:

uBLASFactory.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASFactory

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: uBLASMatrix<Mat>* create_matrix() const
    
        Create empty matrix


    .. cpp:function:: uBLASVector* create_vector() const
    
        Create empty vector


    .. cpp:function:: uBLASVector* create_local_vector() const
    
        Create empty vector (local)


    .. cpp:function:: SparsityPattern* create_pattern() const
    
        Create empty sparsity pattern


    .. cpp:function:: UmfpackLUSolver* create_lu_solver() const
    
        Create LU solver


