
.. Documentation for the header file dolfin/la/PaStiXLUSolver.h

.. _programmers_reference_cpp_la_pastixlusolver:

PaStiXLUSolver.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PaStiXLUSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: PaStiXLUSolver(const STLMatrix& A)
    
        Constructor


    .. cpp:function:: PaStiXLUSolver(std::shared_ptr<const STLMatrix> A)
    
        Constructor


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve Ax = b


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


