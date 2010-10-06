.. Documentation for the header file dolfin/la/EpetraLUSolver.h

.. _programmers_reference_cpp_la_epetralusolver:

EpetraLUSolver.h
================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: EpetraLUSolver

    *Parent class*
    
        * :cpp:class:`GenericLUSolver`
        
    This class implements the direct solution (LU factorization) for
    linear systems of the form Ax = b. It is a wrapper for the LU
    solver of Epetra.

    .. cpp:function:: EpetraLUSolver()
    
        Constructor

    .. cpp:function:: EpetraLUSolver(const GenericMatrix& A)
    
        Constructor

    .. cpp:function:: void set_operator(const GenericMatrix& A)
    
        Set operator (matrix)

    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b

    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b

    .. cpp:function:: uint solve(const EpetraMatrix& A, EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

