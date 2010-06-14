.. Documentation for the Poisson demo from DOLFIN.

.. _demos_python_pde_poisson:

Poisson's equation
==================

.. include:: ../../../common/pde/poisson/poisson.txt


Implementation
--------------

.. code-block:: python

    # Define variational problem
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
    h = Expression("-sin(5*x[0])")
    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*g*ds

