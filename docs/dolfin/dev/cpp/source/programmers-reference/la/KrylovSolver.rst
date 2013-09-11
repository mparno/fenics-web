
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
        
    This class defines an interface for a Krylov solver. The
    approproiate solver is chosen on the basis of the matrix/vector
    type.


    .. cpp:function:: KrylovSolver(std::string method="default", std::string preconditioner="default")
    
        Constructor


    .. cpp:function:: KrylovSolver(boost::shared_ptr<const GenericLinearOperator> A, std::string method="default", std::string preconditioner="default")
    
        Constructor


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericLinearOperator> A, const boost::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& nullspace)
    
        Set null space of the operator (matrix). This is used to solve
        singular systems


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


