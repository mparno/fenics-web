.. Documentation for the header file dolfin/la/ITLKrylovSolver.h

.. _programmers_reference_cpp_la_Mesh:

ITLKrylovSolver.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: class GenericMatrix
    
        Forward declarations

.. cpp:class:: ITLKrylovSolver

    *Parent class*
    
        * :cpp:class:`GenericLinearSolver`
        
        This class implements Krylov methods for linear systems
        of the form Ax = b. It is a wrapper for the Krylov solvers
        of ITL.

    .. cpp:function:: ITLKrylovSolver(std::string method = "default", std::string pc_type = "default")
    
        Create Krylov solver for a particular method and preconditioner

    .. cpp:function:: boost::shared_ptr<const MTL4Matrix> A
    
        Operator (the matrix)

    .. cpp:function:: boost::shared_ptr<const MTL4Matrix> P
    
        Matrix used to construct the preconditoner

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations

    .. cpp:function:: uint solve(MTL4Vector& x, const MTL4Vector& b)
    
        Solve linear system Ax = b and return number of iterations

    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations

    .. cpp:function:: void set_operator(const GenericMatrix& A)
    
        Set operator (matrix)

    .. cpp:function:: void set_operators(const GenericMatrix& A, const GenericMatrix& P)
    
        Set operator (matrix) and preconditioner matrix

    .. cpp:function:: ~ITLKrylovSolver()
    
        Destructor

