
.. Documentation for the header file dolfin/adaptivity/adaptivesolve.h

.. _programmers_reference_cpp_adaptivity_adaptivesolve:

adaptivesolve.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void solve(const Equation& equation, Function& u, const double tol, GoalFunctional& M)

    Solve linear variational problem a(u, v) == L(v) without
    essential boundary conditions


.. cpp:function:: void solve(const Equation& equation, Function& u, const DirichletBC& bc, const double tol, GoalFunctional& M)

    Solve linear variational problem a(u, v) == L(v) with single
    boundary condition


.. cpp:function:: void solve(const Equation& equation, Function& u, std::vector<const DirichletBC*> bcs, const double tol, GoalFunctional& M)

    Solve linear variational problem a(u, v) == L(v) with list of
    boundary conditions


.. cpp:function:: void solve(const Equation& equation, Function& u, const Form& J, const double tol, GoalFunctional& M)

    Solve nonlinear variational problem F(u; v) = 0 without
    essential boundary conditions


.. cpp:function:: void solve(const Equation& equation, Function& u, const DirichletBC& bc, const Form& J, const double tol, GoalFunctional& M)

    Solve linear variational problem F(u; v) = 0 with single
    boundary condition


.. cpp:function:: void solve(const Equation& equation, Function& u, std::vector<const DirichletBC*> bcs, const Form& J, const double tol, GoalFunctional& M)

    Solve linear variational problem F(u; v) = 0 with list of
    boundary conditions


