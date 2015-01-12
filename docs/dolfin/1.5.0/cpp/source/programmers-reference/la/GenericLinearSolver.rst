
.. Documentation for the header file dolfin/la/GenericLinearSolver.h

.. _programmers_reference_cpp_la_genericlinearsolver:

GenericLinearSolver.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericLinearSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class provides a general solver for linear systems Ax = b.


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A) = 0
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(std::shared_ptr<const GenericLinearOperator> A, std::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: void set_nullspace(const VectorSpaceBasis& nullspace)
    
        Set null space of the operator (matrix). This is used to solve
        singular systems


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve_transpose(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve_transpose(GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::string parameter_type() const
    
        Return parameter type: "krylov_solver" or "lu_solver"


    .. cpp:function:: void update_parameters(const Parameters& parameters)
    
        Update solver parameters (useful for LinearSolver wrapper)


