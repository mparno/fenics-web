
.. Documentation for the header file dolfin/la/MUMPSLUSolver.h

.. _programmers_reference_cpp_la_mumpslusolver:

MUMPSLUSolver.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MUMPSLUSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class implements the direct solution (LU factorization) for
    linear systems of the form Ax = b. It is a wrapper for the MUMPS
    LU solver.


    .. cpp:function:: MUMPSLUSolver(const CoordinateMatrix& A)
    
        Constructor


    .. cpp:function:: MUMPSLUSolver(std::shared_ptr<const CoordinateMatrix> A)
    
        Constructor


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


