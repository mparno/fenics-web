
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


    .. cpp:function:: NonlinearVariationalSolver(boost::shared_ptr<NonlinearVariationalProblem> problem)
    
        Create nonlinear variational solver for given problem (shared pointer version)


    .. cpp:function:: void solve()
    
        Solve variational problem


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


.. cpp:class:: NonlinearDiscreteProblem

    *Parent class(es)*
    
        * :cpp:class:`NonlinearProblem`
        
