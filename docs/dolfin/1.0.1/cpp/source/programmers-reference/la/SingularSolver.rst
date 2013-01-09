
.. Documentation for the header file dolfin/la/SingularSolver.h

.. _programmers_reference_cpp_la_singularsolver:

SingularSolver.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SingularSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class provides a linear solver for singular linear systems
    Ax = b where A has a one-dimensional null-space (kernel). This
    may happen for example when solving Poisson's equation with
    pure Neumann boundary conditions.
    
    The solver attempts to create an extended non-singular system
    by adding the constraint [1, 1, 1, ...]^T x = 0.
    
    If an optional mass matrix M is supplied, the solver attempts
    to create an extended non-singular system by adding the
    constraint m^T x = 0 where m is the lumped mass matrix. This
    corresponds to setting the average (integral) of the finite
    element function with coefficients x to zero.
    
    The solver makes not attempt to check that the null-space is
    indeed one-dimensional. It is also assumed that the system
    Ax = b retains its sparsity pattern between calls to solve().


    .. cpp:function:: SingularSolver(std::string method = "lu", std::string preconditioner = "ilu")
    
        Create linear solver


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b, const GenericMatrix& M)
    
        Solve linear system Ax = b using mass matrix M for setting constraint


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


