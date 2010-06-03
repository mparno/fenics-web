.. Documentation for the Poisson demo from DOLFIN.

.. demos_python_pde_poisson:

Poisson's equation
==================

First some common introduction

.. include:: ../../../common/pde/poisson/poisson.txt

This text is Python specific.
Here is a code snippet from the Poisson :download:`demo.py`:

.. code-block:: python

    # Define variational problem
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
    g = Expression("-sin(5*x[0])")
    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*g*ds

