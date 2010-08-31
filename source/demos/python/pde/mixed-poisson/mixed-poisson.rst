.. Documentation for the mixed Poisson demo from DOLFIN.

.. _demos_python_pde_mixed-poisson:

Mixed formulation for Poisson's equation
==========================================

.. include:: ../../../common/pde/mixed-poisson/mixed-poisson.txt


Implementation
--------------

This demo is implemented in a single Python file, :download:`demo.py`, which
contains both the variational forms and the solver.

First, the ``dolfin`` module is imported:

.. code-block:: python

    from dolfin import *

.. index:: FunctionSpace

Then, we need to create a ``mesh`` covering the unit square. In this
example, we will let the mesh consist of 32 x 32 squares with each
square divided into two triangles:

.. code-block:: python

    # Create mesh
    mesh = UnitSquare(32, 32)

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
order. The supported element families are "CG" (continuous elements),
"BDM" (Brezzi-Douglas-Marini), "RT" (Raviart-Thomas), "NED" (Nedelec
elements of the first kind) and "DG" (discontinuous), see the UFL user
manual for more details.  The * operator creates a mixed (product)
space W from the two separate spaces.

Next, we need to specify the trial functions (the unknowns) and the
test functions on this space. This can be done as follows

.. code-block:: python

    # Define trial and test functions
    (sigma, u) = TrialFunctions(V)
    (tau, v) = TestFunctions(V)

In order to define the variational form, it only remains to define
the source function. This is done just as for the Poisson demo:

.. code-block:: python

    # Define source function
    f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")

We are now ready to define the variational forms a and L. Since,
:math:`u_0 = 0` in this example, the boundary term on the right-hand
side vanishes.

.. code-block:: python

    # Define variational form
    a = (dot(sigma, tau) - u*div(tau) + div(sigma)*v)*dx
    L = v*f*dx

It only remains to prescribe the boundary condition for the flux. Need
to construct a :math:`G` such that :math:`G \cdot n = g`. A simple
``Python`` function, which returns a ``bool``, is used to define the
subdomain for the essential boundary condition.  For instance, we can
let

.. code-block:: python

    # Define essential boundary
    def boundary(x):
        return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS


.. code-block:: python

    # Define function G such that G \cdot n = g
    class Flux(Expression):
        def eval_data(self, values, data):
            g = - sin(5*data.x()[0])
            values[0] = g*data.normal()[0]
            values[1] = g*data.normal()[1]
        def rank(self):
            return 1
        def dim(self):
            return 2

    G = Flux()


A Dirichlet boundary condition (``DirichletBC``) can be created. The value
of the boundary condition is represented using a ``Constant``
(equal to zero for the considered case) ), then the  boundary condition is
created:

.. code-block:: python

    bc = DirichletBC(V.sub(0), G, on_boundary)

To compute the solution, a ``VariationalProblem`` object is created
using the bilinear and linear forms, and the boundary condition.  The
``solve`` function is then called, yielding the solution. The separate
components ``sigma`` and ``u`` of the solution can be extracted by
calling the ``split`` function.

.. code-block:: python

    # Compute solution
    problem = VariationalProblem(a, L, bc)
    (sigma, u) = problem.solve().split()

    plot(sigma)
    plot(u, interactive=True)

Complete code
-------------
.. literalinclude:: demo.py
   :start-after: # Begin demo
