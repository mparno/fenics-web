
.. Documentation for the header file dolfin/la/LUSolver.h

.. _programmers_reference_cpp_la_lusolver:

LUSolver.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LUSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLUSolver`
        
    LU solver for the built-in LA backends.


    .. cpp:function:: LUSolver(std::string method= "default")
    
        Constructor


    .. cpp:function:: LUSolver(boost::shared_ptr<const GenericLinearOperator> A, std::string method="default")
    
        Constructor


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve_transpose(GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system


    .. cpp:function:: std::size_t solve_transpose(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void update_parameters(const Parameters& parameters)
    
        Update solver parameters (pass parameters down to wrapped implementation)


