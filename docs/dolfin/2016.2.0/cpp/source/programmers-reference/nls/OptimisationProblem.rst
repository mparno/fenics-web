
.. Documentation for the header file dolfin/nls/OptimisationProblem.h

.. _programmers_reference_cpp_nls_optimisationproblem:

OptimisationProblem.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: OptimisationProblem

    *Parent class(es)*
    
        * :cpp:class:`NonlinearProblem`
        
    This is a base class for nonlinear optimisation problems which
    return the real-valued objective function :math:`f(x)`, its
    gradient :math:`F(x) = f'(x)` and its Hessian :math:`J(x) =
    f''(x)`


    .. cpp:function:: OptimisationProblem()
    
        Constructor


    .. cpp:function:: double f(const GenericVector& x) = 0
    
        Compute the objective function :math:`f(x)`


    .. cpp:function:: void form(GenericMatrix& A, GenericVector& b, const GenericVector& x)
    
        Compute the Hessian :math:`J(x)=f''(x)` and the gradient
        :math:`F(x)=f'(x)`


    .. cpp:function:: void F(GenericVector& b, const GenericVector& x) = 0
    
        Compute the gradient :math:`F(x) = f'(x)`


    .. cpp:function:: void J(GenericMatrix& A, const GenericVector& x) = 0
    
        Compute the Hessian :math:`J(x) = f''(x)`


