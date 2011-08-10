
.. Documentation for the header file dolfin/nls/NewtonSolver.h

.. _programmers_reference_cpp_nls_newtonsolver:

NewtonSolver.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: NewtonSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class defines a Newton solver for nonlinear systems of
    equations of the form :math:`F(x) = 0`.


    .. cpp:function:: NewtonSolver(std::string solver_type="lu", std::string pc_type="default")
    
        Create nonlinear solver with default linear solver and default
        linear algebra backend


    .. cpp:function:: NewtonSolver(GenericLinearSolver& solver, LinearAlgebraFactory& factory)
    
        Create nonlinear solver using provided linear solver and linear algebra
        backend determined by factory
        
        *Arguments*
            solver (:cpp:class:`GenericLinearSolver`)
                The linear solver.
            factory (:cpp:class:`LinearAlgebraFactory`)
                The factory.


    .. cpp:function:: std::pair<uint, bool> solve(NonlinearProblem& nonlinear_function, GenericVector& x)
    
        Solve abstract nonlinear problem :math:`F(x) = 0` for given
        :math:`F` and Jacobian :math:`\dfrac{\partial F}{\partial x}`.
        
        *Arguments*
            nonlinear_function (:cpp:class:`NonlinearProblem`)
                The nonlinear problem.
            x (:cpp:class:`GenericVector`)
                The vector.
        
        *Returns*
            std::pair<uint, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: uint iteration() const
    
        Return Newton iteration number
        
        *Returns*
            uint
                The iteration number.


    .. cpp:function:: double residual() const
    
        Return current residual
        
        *Returns*
            double
                Current residual.


    .. cpp:function:: double relative_residual() const
    
        Return current relative residual
        
        *Returns*
            double
              Current relative residual.


    .. cpp:function:: GenericLinearSolver& linear_solver() const
    
        Return the linear solver
        
        *Returns*
            :cpp:class:`GenericLinearSolver`
                The linear solver.


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values
        
        *Returns*
            :cpp:class:`Parameters`
                Parameter values.


    .. cpp:function:: bool converged(const GenericVector& b, const GenericVector& dx, const NonlinearProblem& nonlinear_problem)
    
        Convergence test


