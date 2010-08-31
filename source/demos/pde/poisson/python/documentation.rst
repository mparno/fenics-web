.. Documentation for the Poisson demo from DOLFIN.

.. _demos_python_pde_poisson:

Poisson equation
================

This demo is implemented in a single Python file, :download:`demo.py`,
which contains both the variational forms and the solver.

.. include:: ../common.txt

Implementation
--------------

First, the ``dolfin`` module is imported:

.. code-block:: python

    from dolfin import *

.. index:: FunctionSpace

Then, a ``mesh`` with 32 vertices in each direction is created
and the finite element ``FunctionSpace`` :math:`V`
is created:

.. code-block:: python

    # Create mesh and define function space
    mesh = UnitSquare(32, 32)
    V = FunctionSpace(mesh, "Lagrange", 1)

The space :math:`V` involves first-order, continuous Lagrange finite element
basis functions.

.. index:: Subdomain

A simple ``Python`` function, which returns a ``bool``, is used to define the
subdomain for the Dirichlet boundary condition:

.. code-block:: python

    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

.. index:: DirichletBC

A Dirichlet boundary condition (``DirichletBC``) can be created. The value
of the boundary condition is represented using a ``Constant``
(equal to zero for the considered case) ), then the  boundary condition is
created:

.. code-block:: python

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)

The space :math:`V` is present to indicate that the Dirichlet boundary
is applied to functions in :math:`V`.

.. index:: Expression

Next, the variational problem is expressed in UFL. First,
``TestFunction`` and ``TrialFunction`` are declared on the already defined
``FunctionSpace`` :math:`V`. Then, the source term :math:`f` and the
normal derivative term on the boundary  :math:`g` are declared
using ``Expression`` class.
Note that the strings defining ``f`` and ``g`` use C++ syntax
since, for efficiency, DOLFIN will generate and compile C++ code for these
expressions at runtime.
The bilinear form ``a`` and the linear form ``L`` are expressed
in UFL.

.. code-block:: python

    # Define variational problem
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
    g = Expression("sin(5*x[0])")
    a = inner(grad(v), grad(u))*dx
    L = v*f*dx + v*g*ds

.. index:: VariationalProblem

To compute the solution, a ``VariationalProblem`` object is created
using the bilinear and linear forms, and the Dirichlet boundary condition.
The ``solve`` function is then called and the solution is returned
in variable ``u`` (a ``Function`` object):.

.. code-block:: python

    # Compute solution
    problem = VariationalProblem(a, L, bc)
    u = problem.solve()

The default settings for solving a variational problem have been used.
Various parameters can be set control aspects of the solution process.

Finally, the solution is output a ``VTK`` file for later visualization
and the the solution plotted to the screen using the ``plot()`` command:

.. code-block:: python

    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    # Plot solution
    plot(u, interactive=True)

Complete code
-------------

.. literalinclude:: demo.py
   :start-after: # Begin demo
