
.. Documentation for the header file dolfin/nls/PETScTAOSolver.h

.. _programmers_reference_cpp_nls_petsctaosolver:

PETScTAOSolver.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScTAOSolver

    *Parent class(es)*
    
        * :cpp:class:`PETScObject`
        
    This class implements methods for solving nonlinear optimisation
    problems via PETSc TAO solver. It supports unconstrained as well
    as bound-constrained minimisation problem


    .. cpp:function:: PETScTAOSolver(const std::string tao_type="default", const std::string ksp_type="default", const std::string pc_type="default")
    
        Create TAO solver for a particular method


    .. cpp:function:: std::size_t solve(OptimisationProblem& optimisation_problem, GenericVector& x, const GenericVector& lb, const GenericVector& ub)
    
        Solve a nonlinear bound-constrained optimisation problem
        
        *Arguments*
            optimisation_problem (:cpp:class:`OptimisationProblem`)
                The nonlinear optimisation problem.
            x (:cpp:class:`GenericVector`)
                The solution vector (initial guess).
            lb (:cpp:class:`GenericVector`)
                The lower bound.
            ub (:cpp:class:`GenericVector`)
                The upper bound.
        
        *Returns*
            num_iterations (std::size_t)
                Number of iterations


    .. cpp:function:: std::size_t solve(OptimisationProblem& optimisation_problem, GenericVector& x)
    
        Solve a nonlinear unconstrained minimisation problem
        
        *Arguments*
            optimisation_problem (:cpp:class:`OptimisationProblem`)
                The nonlinear optimisation problem.
            x (:cpp:class:`GenericVector`)
                The solution vector (initial guess).
        
        *Returns*
            num_iterations (std::size_t)
                Number of iterations


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: Tao tao() const
    
        Return the TAO pointer


    .. cpp:function:: void init(OptimisationProblem& optimisation_problem, PETScVector& x, const PETScVector& lb, const PETScVector& ub)
    
        Initialise the TAO solver for a bound-constrained minimisation
        problem, in case the user wants to access the TAO object
        directly


    .. cpp:function:: void init(OptimisationProblem& optimisation_problem, PETScVector& x)
    
        Initialise the TAO solver for an unconstrained minimisation
        problem, in case the user wants to access the TAO object
        directly


    .. cpp:function:: std::size_t solve(OptimisationProblem& optimisation_problem, PETScVector& x, const PETScVector& lb, const PETScVector& ub)
    
        Solve a nonlinear bound-constrained minimisation problem
        
        *Arguments*
            optimisation_problem (:cpp:class:`OptimisationProblem`)
                The nonlinear optimisation problem.
            x (:cpp:class:`PETScVector`)
                The solution vector (initial guess).
            lb (:cpp:class:`PETScVector`)
                The lower bound.
            ub (:cpp:class:`PETScVector`)
                The upper bound.
        
        *Returns*
            num_iterations (std::size_t)
                Number of iterations


