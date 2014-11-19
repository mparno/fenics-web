
.. Documentation for the header file dolfin/la/uBLASPreconditioner.h

.. _programmers_reference_cpp_la_ublaspreconditioner:

uBLASPreconditioner.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASPreconditioner

    This class specifies the interface for preconditioners for the
    uBLAS Krylov solver.


    .. cpp:function:: uBLASPreconditioner()
    
        Constructor


    .. cpp:function:: void init(const uBLASMatrix<ublas_sparse_matrix>& P)
    
        Initialise preconditioner (sparse matrix)


    .. cpp:function:: void init(const uBLASMatrix<ublas_dense_matrix>& P)
    
        Initialise preconditioner (dense matrix)


    .. cpp:function:: void init(const uBLASLinearOperator& P)
    
        Initialise preconditioner (virtual matrix)


    .. cpp:function:: void solve(uBLASVector& x, const uBLASVector& b) const = 0
    
        Solve linear system (M^-1)Ax = y


