.. Documentation for the Cahn-Hilliard demo from DOLFIN.

.. _demos_cpp_pde_cahn_hilliard:

Cahn-Hilliard equation
======================

.. include:: ../common.txt

This text is C++ specific

UFL form file
-------------

The UFL code for this problem in two dimensions is
:download:`CahnHilliard2D.ufl`.

First, a mixed function spaces of linear Lagrange functions on triangles
is created:

.. code-block:: python

    P1 = FiniteElement("Lagrange", "triangle", 1)
    ME = P1*P1

On the mixed space, test and trial functions are defined:

.. code-block:: python

    q, v = TestFunctions(ME)
    du   = TrialFunction(ME)

The test functions have been split into components.

Coefficient functions are now defined for the current solution (the most recent
guess) and the solution from the beginning of the time step, and the functions
and split into components:

.. code-block:: python

    u   = Coefficient(ME)  # current solution
    u0  = Coefficient(ME)  # solution from previous converged step

    # Split mixed functions
    dc, dmu = split(du)
    c,  mu  = split(u)
    c0, mu0 = split(u0)

Various model parameters are defined as ``Constants``. This means that their
value can be changed without recompiling the UFL file.
Lastly, the value of :math:`\mu_{n+\theta}` is computed.

.. code-block:: python

    lmbda    = Constant(triangle) # surface energy parameter
    dt       = Constant(triangle) # time step
    theta    = Constant(triangle) # time stepping parameter

    # mu_(n+theta)
    mu_mid = (1-theta)*mu0 + theta*mu

To compute :math:`df/dc`, ``c`` is made a variable which will permit
differentiation with respect to it. The function :math:`f` is defined as
a function of :math:`c`, and then

.. code-block:: python

    # Compute the chemical potential df/dc
    c    = variable(c)
    f    = 100*c**2*(1-c)**2
    dfdc = diff(f, c)

The fully discrete variational problem is input, and the Jacobian of the
functional ``L`` which we we wish to drive to zero during the solution process
is computed using the directional derivative (``a``).

.. code-block:: python

    L0 = q*c*dx  - q*c0*dx   + dt*dot(grad(q), grad(mu_mid))*dx
    L1 = v*mu*dx - v*dfdc*dx - lmbda*dot(grad(v), grad(c))*dx
    L = L0 + L1

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

