.. Documentation for the Cahn-Hilliard demo from DOLFIN.

.. _demos_cpp_pde_cahn_hilliard:

Cahn-Hilliard equation
======================

.. include:: ../../../common/pde/cahn-hilliard/cahn-hilliard.txt

This text is C++ specific


UFL input
---------

Some UFL code in :download:`CahnHilliard2D.ufl` and
:download:`CahnHilliard3D.ufl`.

.. code-block:: python

    P1 = FiniteElement("Lagrange", "triangle", 1)
    ME = P1 + P1

    q, v  = TestFunctions(ME)
    du    = TrialFunction(ME)

    u   = Coefficient(ME)  # current solution
    u0  = Coefficient(ME)  # solution from previous converged step

    # Split mixed functions
    dk, dc = split(du)
    k,  c  = split(u)
    k0, c0 = split(u0)

    lmbda    = Constant(triangle) # surface parameter
    muFactor = Constant(triangle) # chemical free energy multiplier

    dt       = Constant(triangle) # time step
    theta    = Constant(triangle) # time stepping parameter

    # Potential mu = \phi,c (chemical free-energy \phi = c^2*(1-c)^2)
    mu = muFactor*(2*c*(1-c)*(1-c) - 2*c*c*(1-c))

    # k^(n+theta)
    k_mid = (1-theta)*k0 + theta*k

    L1 = q*c*dx - q*c0*dx + dt*dot(grad(q), grad(k_mid))*dx
    L2 = v*k*dx - v*mu*dx - lmbda*dot(grad(v), grad(c))*dx
    L = L1 + L2

    a = derivative(L, u, du)


C++ code
--------

.. code-block:: c++

    #include <dolfin.h>
    #include "CahnHilliard2D.h"
    #include "CahnHilliard3D.h"

    using namespace dolfin;

Some more C++

.. code-block:: c++

    // Initial conditions
    class InitialConditions : public Expression
    {
    public:

      InitialConditions(const Mesh& mesh) : Expression(mesh.topology().dim())
      {
        dolfin::seed(2);
      }

      void eval(Array<double>& values, const Data& data) const
      {
        values[0]= 0.0;
        values[1]= 0.63 + 0.02*(0.5 - dolfin::rand());
      }

    };

Yet more C++

.. code-block:: c++

    // Solve
    while (t < T)
    {
      // Update for next time step
      t += dt;
      u0.vector() = u.vector();

      // Solve
      newton_solver.solve(cahn_hilliard, u.vector());

      // Save function to file
      file << u[1];
    }

All this should be in the same :download:`main.cpp` file.

