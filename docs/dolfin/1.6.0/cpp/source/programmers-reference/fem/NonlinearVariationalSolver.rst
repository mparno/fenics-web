
.. Documentation for the header file dolfin/fem/NonlinearVariationalSolver.h

.. _programmers_reference_cpp_fem_nonlinearvariationalsolver:

NonlinearVariationalSolver.h
============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: NonlinearVariationalSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class implements a solver for nonlinear variational problems.


    .. cpp:function:: NonlinearVariationalSolver(NonlinearVariationalProblem& problem)
    
        Create nonlinear variational solver for given problem


    .. cpp:function:: NonlinearVariationalSolver(std::shared_ptr<NonlinearVariationalProblem> problem)
    
        Create nonlinear variational solver for given problem (shared
        pointer version)


    .. cpp:function:: std::pair<std::size_t, bool> solve(const GenericVector& lb, const GenericVector& ub)
    
        Solve variational problem with bound constraints defined by
        GenericVectors
        
        *Arguments*
            lb (:cpp:class:`GenericVector`)
                The linear solver.
            ub (:cpp:class:`GenericVector`)
                The factory.
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::pair<std::size_t, bool> solve(std::shared_ptr<const GenericVector> lb, std::shared_ptr<const GenericVector> ub)
    
        Solve variational problem with bound constraints defined by
        GenericVectors (shared pointer version)
        
        *Arguments*
            lb (_std::shared_ptr<const GenericVector>_)
                The linear solver.
            ub (_std::shared_ptr<const GenericVector>_)
                The factory.
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::pair<std::size_t, bool> solve(const Function& lb, const Function& ub)
    
        Solve variational problem with bound constraints defined by Functions
        
        *Arguments*
            lb (:cpp:class:`Function`)
                The linear solver.
            ub (:cpp:class:`Function`)
                The factory.
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::pair<std::size_t, bool> solve(std::shared_ptr<const Function> lb, std::shared_ptr<const Function> ub)
    
        Solve variational problem with bound constraints defined by
        Functions (shared pointer version)
        
        *Arguments*
            lb (_std::shared_ptr<const Function>_)
                The linear solver.
            ub (_std::shared_ptr<const Function>_)
                The factory.
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::pair<std::size_t, bool> solve()
    
        Solve variational problem
        
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


.. cpp:class:: NonlinearDiscreteProblem

    *Parent class(es)*
    
        * :cpp:class:`NonlinearProblem`
        
