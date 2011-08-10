
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
        
    .. cpp:function:: LUSolver(const GenericMatrix& A, std::string type = "lu")
    
        Constructor


    .. cpp:function:: void set_operator(const GenericMatrix& A)
    
        Set operator (matrix)


    .. cpp:function:: const GenericMatrix& get_operator() const
    
        Return the operator (matrix)


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


