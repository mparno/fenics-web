.. Documentation for the Poisson demo from DOLFIN.

.. _demos_cpp_pde_poisson:


Poisson equation
================

.. include:: ../common.txt


Implementation
--------------

The implementation is split in two files, a form file containing the definition
of the variational forms expressed in UFL and the solver which is implemented
in a C++ file.

Creating the form file
^^^^^^^^^^^^^^^^^^^^^^

First we define the variational problem in UFL which we save in the file called
:download:`Poisson.ufl`.
We first define the finite element, in this case a linear Lagrange triangle:

.. code-block:: python

    element = FiniteElement("Lagrange", triangle, 1)

Then we use this element to initialise the test and trial functions (:math:`v`
and :math:`u`) and coefficient functions (:math:`f` and :math:`g`):

.. code-block:: python

    v = TestFunction(element)
    u = TrialFunction(element)
    f = Coefficient(element)
    g = Coefficient(element)

Finally, we define the bilinear and linear forms according to the equations:

.. code-block:: python

    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*g*ds


Writing the solver
^^^^^^^^^^^^^^^^^^

The solver is implemented in the :download:`main.cpp` file.

At the top we include the DOLFIN header file and the header file containing the
variational forms for the Poisson equation.
For convenience we also include the DOLFIN namespace.

.. code-block:: c++

    #include <dolfin.h>
    #include "Poisson.h"

    using namespace dolfin;

.. index:: Expression

Then follows the definition of the coefficient functions (for :math:`f` and
:math:`g`), which are derived from the ``Expression`` class in DOLFIN.

.. code-block:: c++

    // Source term (right-hand side)
    class Source : public Expression
    {
      void eval(Array<double>& values, const Array<double>& x) const
      {
        double dx = x[0] - 0.5;
        double dy = x[1] - 0.5;
        values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
      }
    };

    // Normal derivative (Neumann boundary condition)
    class dUdN : public Expression
    {
      void eval(Array<double>& values, const Array<double>& x) const
      {
        values[0] = sin(5*x[0]);
      }
    };

.. index:: Subdomain

The ``DirichletBoundary`` is derived from the ``Subdomain`` class and defines
the part of the boundary to which the Dirichlet boundary condition should be
applied.

.. code-block:: c++

    // Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public SubDomain
    {
      bool inside(const Array<double>& x, bool on_boundary) const
      {
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
      }
    };

Inside the ``main()`` function we first create the ``mesh`` and then we define
the ``FunctionSpace`` :math:`V` for our finite element functions.

.. code-block:: c++

    // Create mesh and function space
    UnitSquare mesh(32, 32);
    Poisson::FunctionSpace V(mesh);

.. index:: DirichletBC

After creating the ``FunctionSpace`` and defining our ``DirichletBoundary``
class, we can create the Dirichlet boundary condition (``DirichletBC``) for our
variational problem where we use a ``Constant`` (equal to zero) for the value
of :math:`u` on the Dirichlet boundary.

.. code-block:: c++

    // Define boundary condition
    Constant u0(0.0);
    DirichletBoundary boundary;
    DirichletBC bc(V, u0, boundary);

.. index::
    triple: forms; attach; expression

Next, we define the variational problem by initialising the bilinear and linear
forms (:math:`a`, :math:`L`) using the previously defined ``FunctionSpace``
:math:`V`.
Then we can create the source and boundary flux term (:math:`f`, :math:`g`) and
attach these to the linear form.

.. code-block:: c++

    // Define variational problem
    Poisson::BilinearForm a(V, V);
    Poisson::LinearForm L(V);
    Source f;
    dUdN g;
    L.f = f;
    L.g = g;

.. index:: VariationalProblem

To compute the solution we use the ``VariationalProblem`` class and choose an
iterative linear solver.
The solution is stored in the ``Function`` ``u`` which we also initialise using
the ``FunctionSpace`` :math:`V` since our objective is to find :math:`u \in V`.

.. code-block:: c++

    // Compute solution
    VariationalProblem problem(a, L, bc);
    problem.parameters["linear_solver"] = "iterative";
    Function u(V);
    problem.solve(u);

Finally, we can write the solution to a ``VTK`` file and visualise the solution
using the ``plot()`` command.

.. code-block:: c++

    // Save solution in VTK format
    File file("poisson.pvd");
    file << u;

    // Plot solution
    plot(u);

