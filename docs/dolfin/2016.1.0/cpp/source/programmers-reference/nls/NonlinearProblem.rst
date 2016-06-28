
.. Documentation for the header file dolfin/nls/NonlinearProblem.h

.. _programmers_reference_cpp_nls_nonlinearproblem:

NonlinearProblem.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: NonlinearProblem

    This is a base class for nonlinear problems which can return the
    nonlinear function F(u) and its Jacobian J = dF(u)/du.


    .. cpp:function:: NonlinearProblem()
    
        Constructor


    .. cpp:function:: void form(GenericMatrix& A, GenericVector& b, const GenericVector& x)
    
        Function called by Newton solver before requesting F or J.
        This can be used to compute F and J together


    .. cpp:function:: void F(GenericVector& b, const GenericVector& x) = 0
    
        Compute F at current point x


    .. cpp:function:: void J(GenericMatrix& A, const GenericVector& x) = 0
    
        Compute J = F' at current point x


