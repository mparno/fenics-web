
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
        
    This class implements a solver for nonlinear variational
    problems.


    .. cpp:function:: explicit NonlinearVariationalSolver(std::shared_ptr<NonlinearVariationalProblem> problem)
    
        Create nonlinear variational solver for given problem


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
        
