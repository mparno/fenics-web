.. Documentation for the Poisson demo from DOLFIN.

.. _demos_cpp_pde_poisson:


Poisson's equation
==================

First some common introduction

.. include:: ../../../common/pde/poisson/poisson.txt

This text is C++ specific

Some UFL code in :download:`Poisson.ufl`.

.. code-block:: python

    element = FiniteElement("Lagrange", triangle, 1)

    v = TestFunction(element)
    u = TrialFunction(element)
    f = Coefficient(element)
    g = Coefficient(element)

    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*g*ds


Some C++

.. code-block:: c++

    #include <dolfin.h>
    #include "Poisson.h"

    using namespace dolfin;

Some more C++

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

Yet more C++

.. code-block:: c++

    // Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public SubDomain
    {
      bool inside(const Array<double>& x, bool on_boundary) const
      {
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
      }
    };

All this should be in the same :download:`main.cpp` file.

