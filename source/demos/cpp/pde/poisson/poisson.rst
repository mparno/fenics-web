.. Documentation for the Poisson demo from DOLFIN.

.. _demos_cpp_pde_poisson:


Poisson's equation
==================

.. include:: ../../../common/pde/poisson/poisson.txt


Implementation
--------------

Creating the form file
^^^^^^^^^^^^^^^^^^^^^^

First we define the variational problem in UFL which we save in the file called
:download:`Poisson.ufl`.
We first define the finite element, in this case a linear Lagrange triangle:

.. code-block:: python

    element = FiniteElement("Lagrange", triangle, 1)

Then we use this element to initialise the test and trial functions (:math:`v`
and :math:`u`) and coefficient functions (:math:`f` and :math:`h`):

.. code-block:: python

    v = TestFunction(element)
    u = TrialFunction(element)
    f = Coefficient(element)
    h = Coefficient(element)

Finally, we define the bilinear and linear forms according to the equations:

.. code-block:: python

    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*h*ds


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

Then follows the definition of the coefficient functions (for :math:`f` and
:math:`h`), which are derived from the ``Expression`` class in DOLFIN.

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

    // Boundary flux (Neumann boundary condition)
    class Flux : public Expression
    {
      void eval(Array<double>& values, const Array<double>& x) const
      {
        values[0] = -sin(5*x[0]);
      }
    };

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

Inside the ``main()`` function


