
.. Documentation for the header file dolfin/fem/solve.h

.. _programmers_reference_cpp_fem_solve:

solve.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void solve(const Equation& equation, Function& u)

    Solve linear variational problem a(u, v) == L(v) or nonlinear
    variational problem F(u; v) = 0 without boundary conditions


.. cpp:function:: void solve(const Equation& equation, Function& u, const BoundaryCondition& bc)

    Solve linear variational problem a(u, v) == L(v) or nonlinear
    variational problem F(u; v) = 0 with a single boundary condition


.. cpp:function:: void solve(const Equation& equation, Function& u, std::vector<const BoundaryCondition*> bcs)

    Solve linear variational problem a(u, v) == L(v) or nonlinear
    variational problem F(u; v) = 0 with a list of boundary conditions


.. cpp:function:: void solve(const Equation& equation, Function& u, const Form& J)

    Solve nonlinear variational problem F(u; v) == 0 without boundary
    conditions. The argument J should provide the Jacobian bilinear
    form J = dF/du.


.. cpp:function:: void solve(const Equation& equation, Function& u, const BoundaryCondition& bc, const Form& J)

    Solve nonlinear variational problem F(u; v) == 0 with a single
    boundary condition. The argument J should provide the Jacobian
    bilinear form J = dF/du.


.. cpp:function:: void solve(const Equation& equation, Function& u, std::vector<const BoundaryCondition*> bcs, const Form& J)

    Solve nonlinear variational problem F(u; v) == 0 with a list of
    boundary conditions. The argument J should provide the Jacobian
    bilinear form J = dF/du.


