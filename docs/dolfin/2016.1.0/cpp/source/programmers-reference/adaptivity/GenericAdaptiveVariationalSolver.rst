
.. Documentation for the header file dolfin/adaptivity/GenericAdaptiveVariationalSolver.h

.. _programmers_reference_cpp_adaptivity_genericadaptivevariationalsolver:

GenericAdaptiveVariationalSolver.h
==================================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericAdaptiveVariationalSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    An abstract class for goal-oriented adaptive solution of
    variational problems.
    


    .. cpp:function:: void solve(const double tol)
    
        Solve such that the functional error is less than the given
        tolerance. Note that each call to solve is based on the
        leaf-node of the variational problem
        
        *Arguments*
            tol  (double)
                The error tolerance


    .. cpp:function:: std::shared_ptr<const Function> solve_primal() = 0
    
        Solve the primal problem. Must be overloaded in subclass.
        
        *Returns*
            :cpp:class:`Function`
                The solution to the primal problem


    .. cpp:function:: std::vector<std::shared_ptr<const DirichletBC> > extract_bcs() const = 0
    
        Extract the boundary conditions for the primal problem. Must
        be overloaded in subclass.
        
        *Returns*
            std::vector<:cpp:class:`DirichletBC`>
                The primal boundary conditions


    .. cpp:function:: double evaluate_goal(Form& M, std::shared_ptr<const Function> u) const = 0
    
        Evaluate the goal functional. Must be overloaded in subclass.
        
        *Arguments*
           M (:cpp:class:`Form`)
               The functional to be evaluated
           u (:cpp:class:`Function`)
               The function of which to evaluate the functional
        
        *Returns*
            double
                The value of M evaluated at u


    .. cpp:function:: void adapt_problem(std::shared_ptr<const Mesh> mesh) = 0
    
        Adapt the problem to other mesh. Must be overloaded in subclass.
        
        *Arguments*
           mesh (:cpp:class:`Mesh`)
               The other mesh


    .. cpp:function:: std::vector<std::shared_ptr<Parameters> > adaptive_data() const
    
        Return stored adaptive data
        
        *Returns*
           std::vector<:cpp:class:`Parameters`>
               The data stored in the adaptive loop


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values:
        
            "max_iterations"     (int)
            "max_dimension"      (int)
            "plot_mesh"          (bool)
            "save_data"          (bool)
            "data_label"         (std::string)
            "reference"          (double)
            "marking_strategy"   (std::string)
            "marking_fraction"   (double)


    .. cpp:function:: void summary()
    
        Present summary of all adaptive data and parameters


    .. cpp:function:: std::size_t num_dofs_primal() = 0
    
        Return the number of degrees of freedom for primal problem
        
        *Returns*
            _std::size_t_
                The number of degrees of freedom


