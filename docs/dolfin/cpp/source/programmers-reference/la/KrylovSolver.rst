
.. Documentation for the header file dolfin/la/KrylovSolver.h

.. _programmers_reference_cpp_la_krylovsolver:

KrylovSolver.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: KrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class defines an interface for a Krylov solver. The approproiate solver
    is chosen on the basis of the matrix/vector type.


    .. cpp:function:: KrylovSolver(std::string solver_type = "default", std::string pc_type = "default")
    
        Create Krylov solver


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericMatrix> A, const boost::shared_ptr<const GenericMatrix> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


