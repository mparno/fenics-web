
.. Documentation for the header file dolfin/la/LinearSolver.h

.. _programmers_reference_cpp_la_linearsolver:

LinearSolver.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LinearSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class provides a general solver for linear systems Ax = b.


    .. cpp:function:: LinearSolver(std::string method = "default", std::string preconditioner = "default")
    
        Create linear solver


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Set the operator (matrix)


    .. cpp:function:: void set_operators(std::shared_ptr<const GenericLinearOperator> A, std::shared_ptr<const GenericLinearOperator> P)
    
        Set the operator (matrix) and preconditioner matrix


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void update_parameters(const Parameters& parameters)
    
        Update solver parameters (pass parameters down to wrapped implementation)


    .. cpp:function:: std::string parameter_type() const
    
        Return parameter type: "krylov_solver" or "lu_solver"


