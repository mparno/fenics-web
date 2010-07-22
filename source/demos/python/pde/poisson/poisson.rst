.. Documentation for the Poisson demo from DOLFIN.

.. _demos_python_pde_poisson:

Poisson's equation
==================

.. include:: ../../../common/pde/poisson/poisson.txt


Implementation
--------------

The demo is implemented in a single Python file, :download:`demo.py`, which
contains both the variational forms and the solver.

First we import everyting from the ``dolfin`` module.

.. code-block:: python

    from dolfin import * 

.. index:: FunctionSpace

Then we create the ``mesh`` and define the ``FunctionSpace`` :math:`V` for our
finite element functions.

.. code-block:: python

    # Create mesh and define function space
    mesh = UnitSquare(32, 32)
    V = FunctionSpace(mesh, "CG", 1)

.. index:: Subdomain

We use a simple ``Python`` function, which returns a ``bool``, to define the
subdomain for the Dirichlet boundary condition.

.. code-block:: python

    # Define Dirichlet boundary (x = 0 or x = 1)
    def boundary(x):
        return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

.. index:: DirichletBC

We can then create the Dirichlet boundary condition (``DirichletBC``) for our
variational problem where we use a ``Constant`` (equal to zero) for the value
of :math:`u` on the Dirichlet boundary.

.. code-block:: python

    # Define boundary condition
    u0 = Constant(0.0)
    bc = DirichletBC(V, u0, boundary)

.. index:: Expression

Next, we define the variational problem in UFL. First we initialise the
``TestFunction`` and ``TrialFunction`` using the previously defined
``FunctionSpace`` :math:`V`. Then we use the ``Expression`` class to define
the source term and normal derivative term on the boundary
(:math:`f`, :math:`g`).
Note that the string defining these expressions uses ``C++`` syntax.
The bilinear and linear forms (``a`` and ``L``) can then be expressed directly
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

To compute the solution we use the ``VariationalProblem`` class which we
initialise with the bilinear and linear forms and the Dirichlet boundary
conditions.
The solution is stored in variable ``u`` (a ``Function`` object) and we obtain
it by calling the ``solve`` method of the ``VariationalProblem``.

.. code-block:: python

    # Compute solution
    problem = VariationalProblem(a, L, bc)
    u = problem.solve()

Finally, we can write the solution to a ``VTK`` file and visualise the solution
using the ``plot()`` command.

.. code-block:: python

    # Save solution in VTK format
    file = File("poisson.pvd")
    file << u

    # Plot solution
    plot(u, interactive=True)

