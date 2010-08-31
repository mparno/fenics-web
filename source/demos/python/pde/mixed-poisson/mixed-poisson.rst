.. Documentation for the mixed Poisson demo from DOLFIN.

.. _demos_python_pde_mixed-poisson:

Mixed formulation for Poisson's equation
========================================

.. include:: ../../../common/pde/mixed-poisson/mixed-poisson.txt


Implementation
--------------

This demo is implemented in a single Python file, :download:`demo.py`, which
contains both the variational forms and the solver.

First, the ``dolfin`` module is imported:

.. code-block:: python

    from dolfin import *

Then, we need to create a ``mesh`` covering the unit square. In this
example, we will let the mesh consist of 32 x 32 squares with each
square divided into two triangles:

.. code-block:: python

    # Create mesh
    mesh = UnitSquare(32, 32)

.. index::
   pair: FunctionSpace; Brezzi-Douglas-Marini
   pair: FunctionSpace; Discontinous Lagrange

Next, we need to define the function spaces. We define the two
function spaces :math:`\Sigma_h` and :math:`V_h` separately, before
combining these into a mixed function space:

.. code-block:: python

    # Define function spaces and mixed (product) space
    BDM = FunctionSpace(mesh, "BDM", 1)
    DG = FunctionSpace(mesh, "DG", 0)
    W = BDM * DG

The second argument to ``FunctionSpace`` specifies the type of finite
element family, while the third argument specifies the polynomial
degree. The UFL user manual contains a list of all available finite
element families and more details.  The * operator creates a mixed
(product) space ``W`` from the two separate spaces ``BDM`` and
``DG``. Hence,

.. math::

    W = \{ (\tau, v) \ \text{such that} \ \tau \in BDM, v \in DG \}.

Next, we need to specify the trial functions (the unknowns) and the
test functions on this space. This can be done as follows

.. code-block:: python

    # Define trial and test functions
    (sigma, u) = TrialFunctions(W)
    (tau, v) = TestFunctions(W)

In order to define the variational form, it only remains to define the
source function :math:`f`. This is done just as for the Poisson demo:

.. code-block:: python

    # Define source function
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")

We are now ready to define the variational forms a and L. Since,
:math:`u_0 = 0` in this example, the boundary term on the right-hand
side vanishes.

.. code-block:: python

    # Define variational form
    a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
    L = - f*v*dx


It only remains to prescribe the boundary condition for the
flux. Essential boundary conditions are specified through the class
``DirichletBC`` which takes three arguments: the function space the
boundary condition is supposed to be applied to, the data for the
boundary condition, and the relevant part of the boundary.

We want to apply the boundary condition to the first subspace of the
mixed space. This space can be accessed by ``W.sub(0)``. (Do *not* use
the separate space ``BDM`` as this would mess up the numbering.)

Next, we need to construct the data for the boundary condition. An
essential boundary condition is handled by replacing degrees of
freedom by the degrees of freedom evaluated at the given data. The
:math:`BDM` finite element spaces are vector-valued spaces and hence
the degrees of freedom act on vector-valued objects. The effect is
that the user is required to construct a :math:`G` such that :math:`G
\cdot n = g`.  Such a :math:`G` can be constructed by letting :math:`G
= g n`. In particular, it can be easily created by subclassing the
``Expression`` class. Overloading the ``eval_data`` method (instead of
the usual ``eval``) allows us to extract more geometry information
such as the facet normals. Since this is a vector-valued expression,
the methods ``rank`` and ``dim`` must also be specified.

.. index:: Expression

.. code-block:: python

    # Define function G such that G \cdot n = g
    class BoundarySource(Expression):
        def eval_data(self, values, data):
            g = sin(5*data.x()[0])
            values[0] = g*data.normal()[0]
            values[1] = g*data.normal()[1]
        def rank(self):
            return 1
        def dim(self):
            return 2
    G = BoundarySource()

Specifying the relevant part of the boundary can be done as for the
Poisson demo (but now the top and bottom of the unit square is the
essential boundary):

.. code-block:: python

    # Define essential boundary
    def boundary(x):
        return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

Now, all the pieces are in place for the construction of the essential
boundary condition:

.. code-block:: python

    bc = DirichletBC(W.sub(0), G, on_boundary)

To compute the solution, a ``VariationalProblem`` object is created
using the bilinear and linear forms, and the boundary condition.  The
``solve`` function is then called, yielding the full solution. The
separate components ``sigma`` and ``u`` of the solution can be
extracted by calling the ``split`` function. Finally, we plot the
solutions to examine the result.

.. index:: split functions

.. code-block:: python

    # Compute solution
    problem = VariationalProblem(a, L, bc)
    (sigma, u) = problem.solve().split()

    # Plot sigma and u
    plot(sigma)
    plot(u)
    interactive()

Complete code
-------------
.. literalinclude:: demo.py
   :start-after: # Begin demo
