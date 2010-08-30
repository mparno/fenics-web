.. Documentation for the header file dolfin/la/GenericLUSolver.h

.. _programmers_reference_cpp_la_genericlusolver:

GenericLUSolver.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: class GenericVector
    
        Forward declarations

.. cpp:class:: GenericLUSolver

    *Parent class*
    
        * :cpp:class:`GenericLinearSolver`
        
    This a base class for LU solvers

    .. cpp:function:: virtual uint solve(GenericVector& x, const GenericVector& b) = 0
    
        Solve linear system Ax = b

    .. cpp:function:: virtual uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b

    .. cpp:function:: virtual void set_operator(const GenericMatrix& A) = 0
    
        Set operator (matrix)

