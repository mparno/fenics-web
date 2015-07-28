
.. Documentation for the header file dolfin/fem/LinearVariationalSolver.h

.. _programmers_reference_cpp_fem_linearvariationalsolver:

LinearVariationalSolver.h
=========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LinearVariationalSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class implements a solver for linear variational problems.


    .. cpp:function:: LinearVariationalSolver(LinearVariationalProblem& problem)
    
        Create linear variational solver for given problem


    .. cpp:function:: LinearVariationalSolver(std::shared_ptr<LinearVariationalProblem> problem)
    
        Create linear variational solver for given problem (shared
        pointer version)


    .. cpp:function:: void solve()
    
        Solve variational problem


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


