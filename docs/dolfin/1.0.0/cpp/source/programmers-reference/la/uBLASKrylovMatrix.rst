
.. Documentation for the header file dolfin/la/uBLASKrylovMatrix.h

.. _programmers_reference_cpp_la_ublaskrylovmatrix:

uBLASKrylovMatrix.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASKrylovMatrix

    This class provides an interface for matrices that define linear
    systems for the uBLASKrylovSolver. This interface is implemented
    by the classes uBLASSparseMatrix and DenseMatrix. Users may also
    overload the mult() function to specify a linear system only in
    terms of its action.


    .. cpp:function:: uBLASKrylovMatrix()
    
        Constructor


    .. cpp:function:: uint size(uint dim) const = 0
    
        Return number of rows (dim = 0) or columns (dim = 1)


    .. cpp:function:: void mult(const uBLASVector& x, uBLASVector& y) const = 0
    
        Compute product y = Ax


    .. cpp:function:: void solve(uBLASVector& x, const uBLASVector& b)
    
        Solve linear system Ax = b for a Krylov matrix using uBLAS and dense matrices


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


