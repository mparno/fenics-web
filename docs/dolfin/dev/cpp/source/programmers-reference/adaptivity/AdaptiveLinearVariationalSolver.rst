
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


    .. cpp:function:: AdaptiveLinearVariationalSolver(LinearVariationalProblem& problem, GoalFunctional& goal)
    
        Create AdaptiveLinearVariationalSolver
        
        *Arguments*
            problem (:cpp:class:`LinearVariationalProblem`)
                The primal problem
            goal (:cpp:class:`GoalFunctional`)
                The goal functional


    .. cpp:function:: AdaptiveLinearVariationalSolver(boost::shared_ptr<LinearVariationalProblem> problem, boost::shared_ptr<GoalFunctional> goal)
    
        Create AdaptiveLinearVariationalSolver (shared ptr version)
        
        *Arguments*
            problem (:cpp:class:`LinearVariationalProblem`)
                The primal problem
            goal (:cpp:class:`GoalFunctional`)
                The goal functional


    .. cpp:function:: AdaptiveLinearVariationalSolver(boost::shared_ptr<LinearVariationalProblem> problem, boost::shared_ptr<Form> goal, boost::shared_ptr<ErrorControl> control)
    
        Create AdaptiveLinearVariationalSolver from variational
        problem, goal form and error control instance
        
        *Arguments*
            problem (:cpp:class:`LinearVariationalProblem`)
                The primal problem
            goal (:cpp:class:`Form`)
                The goal functional
            control (:cpp:class:`ErrorControl`)
                An error controller object


    .. cpp:function:: boost::shared_ptr<const Function> solve_primal()
    
        Solve the primal problem.
        
        *Returns*
            :cpp:class:`Function`
                The solution to the primal problem


    .. cpp:function:: std::vector<boost::shared_ptr<const DirichletBC> > extract_bcs() const
    
        Extract the boundary conditions for the primal problem.
        
        *Returns*
            std::vector<:cpp:class:`DirichletBC`>
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


    .. cpp:function:: std::size_t num_dofs_primal()
    
        Return the number of degrees of freedom for primal problem
        
        *Returns*
            _std::size_t_
                The number of degrees of freedom


    .. cpp:function:: void init(boost::shared_ptr<LinearVariationalProblem> problem, boost::shared_ptr<GoalFunctional> goal)
    
        Helper function for instance initialization
        
        *Arguments*
           problem (:cpp:class:`LinearVariationalProblem`)
               The primal problem
           u (:cpp:class:`GoalFunctional`)
               The goal functional


