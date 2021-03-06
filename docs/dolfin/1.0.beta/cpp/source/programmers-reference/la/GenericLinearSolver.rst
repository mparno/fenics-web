
.. Documentation for the header file dolfin/la/GenericLinearSolver.h

.. _programmers_reference_cpp_la_genericlinearsolver:

GenericLinearSolver.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericLinearSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class provides a general solver for linear systems Ax = b.


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericMatrix> A) = 0
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericMatrix> A, const boost::shared_ptr<const GenericMatrix> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


