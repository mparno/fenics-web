
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


    .. cpp:function:: NewtonSolver()
    
        Create nonlinear solver


    .. cpp:function:: NewtonSolver(std::shared_ptr<GenericLinearSolver> solver, GenericLinearAlgebraFactory& factory)
    
        Create nonlinear solver using provided linear solver
        
        *Arguments*
            solver (:cpp:class:`GenericLinearSolver`)
                The linear solver.
            factory (:cpp:class:`GenericLinearAlgebraFactory`)
                The factory.


    .. cpp:function:: std::pair<std::size_t, bool> solve(NonlinearProblem& nonlinear_function, GenericVector& x)
    
        Solve abstract nonlinear problem :math:`F(x) = 0` for given
        :math:`F` and Jacobian :math:`\dfrac{\partial F}{\partial x}`.
        
        *Arguments*
            nonlinear_function (:cpp:class:`NonlinearProblem`)
                The nonlinear problem.
            x (:cpp:class:`GenericVector`)
                The vector.
        
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::size_t iteration() const
    
        Return Newton iteration number
        
        *Returns*
            std::size_t
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


