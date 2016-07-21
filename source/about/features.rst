.. _features:

########
Features
########

FEniCS comes packed with features for the computational scientist.
Partial differential equations can be specified in near-**mathematical
notation** (as finite element variational problems) and **solved
automatically**. FEniCS also provides a large library of important
tools for the numerical analyst who wishes to explore and develop new
methods.

******************************************
Automated solution of variational problems
******************************************

Finite element variational problems may be specified in
near-mathematical notation **directly as part of your program**.  For
example, the variational problem for the Poisson equation,

.. math::
   \underbrace{\int_{\Omega} \nabla u \cdot \nabla v \, {\rm d} x}_{a(u, v)}
   =
   \underbrace{\int_{\Omega} f v \, {\rm d} x}_{L(v)}
   \quad \forall v \in V,

can be directly translated to the following FEniCS program:

.. code-block:: python

    u = TrialFunction(V)
    v = TestFunction(V)

    a = dot(grad(u), grad(v))*dx
    L = f*v*dx

Variational problems like the one above may be **solved
automatically** in FEniCS by calling the ``solve()`` function:

.. code-block:: python

    u = Function(V)
    solve(a == L, u, bc)

Automated solution of variational problems is not limited to linear
problems. FEniCS also supports **general nonlinear variational
problems**:

.. math::
   F(u; v) = 0 \quad \forall v \in V.

The solution is automatically computed by Newton's method through
**automatic differentiation**:

.. code-block:: python

    solve(F == 0, u, bc)

**************************************
Automated error control and adaptivity
**************************************

Say you want to solve the above problem adaptively with **automated
control of the error** in the computed solution... No problem, just
specify a *goal functional* :math:`\mathcal{M} : V \rightarrow
\mathbb{R}` (a global scalar functional of your computed solution) and
a tolerance :math:`\epsilon > 0`:

.. code-block:: python

    solve(F == 0, u, bc, tol=epsilon, M=M)

***************************************
An extensive library of finite elements
***************************************

FEniCS provides an extensive library of finite elements. You will find
the standard **Lagrange** elements, but also support for **DG
methods**, vector elements like the **BDM**, **RT** and **Nedelec**
elements, and special element types like the **Crouzeix-Raviart**
element.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <img src="../_images/elements.png" class="img-responsive" alt="Elements">
    </div>
  </div>

*******************************
High performance linear algebra
*******************************

FEniCS provides unified access to **a range of linear algebra
libraries** through a common wrapper layer. Currently supported linear
algebra backends include `PETSc <http://www.mcs.anl.gov/petsc/>`_,
`Trilinos/Tpetra <https://trilinos.org/>`_ and `Eigen3
<http://eigen.tuxfamily.org/>`_.  The backend may be easily switched
by changing the value of a parameter. **Parallel computing** is
supported through the PETSc and Trilinos backends.

********************
Computational meshes
********************

FEniCS provides **fully distributed** simplex meshes in one
(intervals), two (triangles) and three (tetrahedra) space dimensions.
Mesh data may be easily accessed through **mesh iterators**. Meshes
may be **refined adaptively**, and **mesh partitioning** for parallel
computing is supported through interfaces to `SCOTCH
<http://www.labri.fr/perso/pelegrin/scotch/>`_ and `ParMETIS
<http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <img src="../_images/meshes.png" class="img-responsive" alt="Meshes">
    </div>
  </div>

**************
Postprocessing
**************

FEniCS provides built-in plotting for quick and easy inspection of
solutions and meshes. Just call the ``plot()`` command for **live
plotting** of your data:

.. code-block:: python

    plot(mesh)
    plot(u)

You can even plot derived quantities like the gradient of a function:

.. code-block:: python

    plot(grad(u))

For more **advanced postprocessing**, FEniCS provides easy output in
VTK format for visualization in `ParaView <http://www.paraview.org/>`_
or `MayaVi <http://mayavi.sourceforge.net/>`_.


*****************
Language bindings
*****************

FEniCS can be used from both **Python** and **C++**. The two interaces
are very similar and provide the same features (with some small
exceptions). Which interface to choose is a matter of taste, but the
Python interface is easier to work with if you are not already a
seasoned C++ programmer.


***********************
Extensive documentation
***********************

FEniCS comes with **extensive documentation**, consisting of a
:ref:`comprehensive tutorial <tutorial>`, detailed :ref:`API
documentation <documentation>` and a range of :ref:`documented demos
<documentation>`. In addition, the :ref:`700-page FEniCS book <book>`
documents the methodology behind the FEniCS Project and highlights a
number of applications in computational science based on FEniCS.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <img src="../_images/documentation.png" class="img-responsive" alt="Fenics Documentation">
    </div>
  </div>
