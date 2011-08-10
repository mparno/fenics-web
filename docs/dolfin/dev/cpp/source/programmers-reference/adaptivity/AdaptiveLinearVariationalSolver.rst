
.. Documentation for the header file dolfin/adaptivity/AdaptiveLinearVariationalSolver.h

.. _programmers_reference_cpp_adaptivity_adaptivelinearvariationalsolver:

AdaptiveLinearVariationalSolver.h
=================================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: AdaptiveLinearVariationalSolver

    A class for goal-oriented adaptive solution of linear
    variational problems.
    
    For a linear variational problem of the form: find u in V
    satisfying
    
        a(u, v) = L(v) for all v in :math:`\hat V`
    
    and a corresponding conforming discrete problem: find u_h in V_h
    satisfying
    
        a(u_h, v) = L(v) for all v in :math:`\hat V_h`
    
    and a given goal functional M and tolerance tol, the aim is to
    find a V_H and a u_H in V_H satisfying the discrete problem such
    that
    
        \|M(u) - M(u_H)\| < tol
    
    This strategy is based on dual-weighted residual error
    estimators designed and automatically generated for the primal
    problem and subsequent h-adaptivity.


    .. cpp:function:: AdaptiveLinearVariationalSolver(LinearVariationalProblem& problem)
    
        Create AdaptiveLinearVariationalSolver
        
        *Arguments*
            problem (:cpp:class:`LinearVariationalProblem`)
                The primal problem


    .. cpp:function:: AdaptiveLinearVariationalSolver(boost::shared_ptr<LinearVariationalProblem> problem)
    
        Create AdaptiveLinearVariationalSolver
        
        *Arguments*
            problem (:cpp:class:`LinearVariationalProblem`)
                The primal problem


    .. cpp:function:: void solve(const double tol, GoalFunctional& M)
    
        Solve problem such that the error measured in the goal
        functional 'M' is less than the given tolerance using the
        GoalFunctional's ErrorControl object.
        
        *Arguments*
            tol  (double)
                The error tolerance
            goal  (:cpp:class:`GoalFunctional`)
                The goal functional


    .. cpp:function:: boost::shared_ptr<const Function> solve_primal()
    
        Solve the primal problem.
        
        *Returns*
            :cpp:class:`Function`
                The solution to the primal problem


    .. cpp:function:: std::vector<boost::shared_ptr<const BoundaryCondition> > extract_bcs() const
    
        Extract the boundary conditions for the primal problem.
        
        *Returns*
            std::vector<:cpp:class:`BoundaryCondition`>
                The primal boundary conditions


    .. cpp:function:: double evaluate_goal(Form& M, boost::shared_ptr<const Function> u) const
    
        Evaluate the goal functional.
        
        *Arguments*
           M (:cpp:class:`Form`)
               The functional to be evaluated
           u (:cpp:class:`Function`)
               The function at which to evaluate the functional
        
        *Returns*
            double
                The value of M evaluated at u


    .. cpp:function:: void adapt_problem(boost::shared_ptr<const Mesh> mesh)
    
        Adapt the problem to other mesh.
        
        *Arguments*
           mesh (:cpp:class:`Mesh`)
               The other mesh


