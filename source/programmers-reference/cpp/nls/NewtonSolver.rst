.. Documentation for the header file dolfin/nls/NewtonSolver.h

.. _programmers_reference_cpp_nls_newtonsolver:

NewtonSolver.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: NewtonSolver

    *Parent class*
    
        * :cpp:class:`Variable`
        
    This class defines a Newton solver for equations of the form F(u) = 0.

    .. cpp:function:: NewtonSolver(std::string solver_type = "lu", std::string pc_type = "default")
    
        Create nonlinear solver with default linear solver and default
        linear algebra backend

    .. cpp:function:: NewtonSolver(GenericLinearSolver& solver, LinearAlgebraFactory& factory)
    
        Create nonlinear solver using provided linear solver and linear algebra
        backend determined by factory

    .. cpp:function:: std::pair<uint, bool> solve(NonlinearProblem& nonlinear_function, GenericVector& x)
    
        Solve abstract nonlinear problem F(x) = 0 for given vector F and
        Jacobian dF/dx

    .. cpp:function:: uint iteration() const
    
        Return Newton iteration number

    .. cpp:function:: double residual() const
    
        Return current residual

    .. cpp:function:: double relative_residual() const
    
        Return current relative residual

    .. cpp:function:: GenericLinearSolver& linear_solver() const
    
        Return the linear solver

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: bool converged(const GenericVector& b, const GenericVector& dx, const NonlinearProblem& nonlinear_problem)
    
        Convergence test

    .. cpp:function:: uint newton_iteration
    
        Current number of Newton iterations

    .. cpp:function:: double _residual, residual0
    
        Most recent residual and intitial residual

    .. cpp:function:: boost::shared_ptr<GenericLinearSolver> solver
    
        Solver

    .. cpp:function:: boost::scoped_ptr<GenericMatrix> A
    
        Jacobian matrix

    .. cpp:function:: boost::scoped_ptr<GenericVector> dx
    
        Solution vector

    .. cpp:function:: boost::scoped_ptr<GenericVector> b
    
        Resdiual vector

