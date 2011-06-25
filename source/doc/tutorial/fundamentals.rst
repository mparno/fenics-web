.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)

.. _tut:fundamentals:

Fundamentals
============

FEniCS is a user-friendly tool for solving partial differential
equations (PDEs). The goal of this tutorial is get you started with
FEniCS through a series of
simple examples that demonstrate

  * how to define the PDE problem in terms of a variational problem

  * how to define simple domains

  * how to deal with Dirichlet, Neumann, and Robin conditions

  * how to deal with variable coefficients

  * how to deal with domains built of several materials (subdomains)

  * how to compute derived quantities like the flux vector field or
    a functional of the solution

  * how to quickly visualize the mesh, the solution, the flux, etc.

  * how to solve nonlinear PDEs in various ways

  * how to deal with time-dependent PDEs

  * how to set parameters governing solution methods for linear systems

  * how to create domains of more complex shape

The mathematics of the illustrations is kept simple to better focus
on FEniCS functionality and syntax. This means that we mostly use
the Poisson equation and the time-dependent diffusion equation
as model problems, often with input data adjusted such that we get
a very simple solution that can be exactly reproduced by any standard
finite element method over a uniform, structured mesh. This
latter property greatly enhances the verification of the implementations.
Occasionally we insert a physically more relevant example
to remind the reader that changing the PDE and boundary
conditions to something more real might often be a trivial task.

.. With the fundamentals explained, we move on to physically more
.. complicated problems, including systems of PDEs, and show how to build
.. more complete simulation codes.

FEniCS may seem to require a thorough understanding of the abstract
mathematical version of the finite element method as well as
familiarity with the Python programming language.  Nevertheless, it
turns out that many are able to pick up the fundamentals of finite
elements *and* Python programming as they go along with this
tutorial. Simply keep on reading and try out the examples. You will be
amazed of how easy it is to solve PDEs with FEniCS!

Reading this tutorial obviously requires access to a machine where the
FEniCS software is installed. The section :ref:`tut:app:install` explains
briefly how to install the necessary tools.

.. _tut:poisson1:bvp:

The Poisson equation
--------------------

.. index:: Poisson's equation


Our first example regards the Poisson problem,

.. math::

        - \nabla^2 u(\pmb{x}) &= f(\pmb{x}),\quad \pmb{x}\mbox{ in } \Omega,
        \\
        u(\pmb{x}) &= u_0(\pmb{x}),\quad \pmb{x}\mbox{ on } \partial \Omega\thinspace .


Here, :math:`u(\pmb{x})` is the unknown function, :math:`f(\pmb{x})` is a
prescribed function of space, :math:`\nabla^2` is the Laplace operator (also
often written as :math:`\Delta`), :math:`\Omega` is the spatial domain, and
:math:`\partial\Omega` is the boundary of :math:`\Omega`. A stationary PDE like
this, together with a complete set of boundary conditions, constitute
a *boundary-value problem*, which must be precisely stated before
it makes sense to start solving it with FEniCS.

In two space dimensions with coordinates :math:`x` and :math:`y`, we can write out
the Poisson equation in detail:

.. math::


        - {\partial^2 u\over\partial x^2} -
        {\partial^2 u\over\partial y^2} = f(x,y)\thinspace .

The unknown :math:`u` is now a function of two variables, :math:`u(x,y)`, defined
over a two-dimensional domain :math:`\Omega`.

The Poisson equation arises in numerous physical contexts, including
heat conduction, electrostatics, diffusion of substances, twisting of
elastic rods, inviscid fluid flow, and water waves. Moreover, the
equation appears in numerical splitting strategies of more complicated
systems of PDEs, in particular the Navier-Stokes equations.

Solving a physical problem with FEniCS consists
of the following steps:

 1. Identify the PDE and its boundary conditions.

 2. Reformulate the PDE problem as a variational problem.

 3. Make a Python program where the formulas in the variational
    problem are coded, along with definitions of input data such as
    :math:`f`, :math:`u_0`, and a mesh for the spatial domain :math:`\Omega`.

 4. Add statements in the program for solving the variational
    problem, computing derived quantities such as :math:`\nabla u`, and
    visualizing the results.

We shall now go through steps 2--4 in detail.  The key feature of
FEniCS is that steps 3 and 4 result in fairly short code, while most
other software frameworks for PDEs require much more code and more
technically difficult programming.


.. _tut:poisson1:varform:

Variational Formulation
-----------------------

.. index:: variational formulation

FEniCS makes it easy to solve PDEs if finite elements are used for
discretization in space and the problem is expressed as a
*variational problem*. Readers who are not familiar with
variational problems will get a brief introduction to the topic in
this tutorial, but getting and reading
a proper book on the finite element method in addition is encouraged. The section :ref:`tut:appendix:books` contains a list of some suitable
books.

.. index:: test function

.. index:: trial function

The core of the recipe for turning a PDE into a variational problem
is to multiply the PDE by a function :math:`v`, integrate the resulting
equation over :math:`\Omega`, and perform integration by parts of terms with
second-order derivatives. The function :math:`v` which multiplies the PDE
is in the mathematical finite element literature
called a *test function*. The unknown function :math:`u` to be approximated
is referred to
as a *trial function*. The terms test and trial function are used
in FEniCS programs too.
Suitable
function spaces must be specified for the test and trial functions.
For standard PDEs arising in physics and mechanics such spaces are
well known.

In the present case, we first multiply the Poisson equation
by the test function :math:`v` and integrate,

.. math::

         -\int_\Omega (\Delta u)v dx = \int_\Omega fv dx\thinspace .

Then we apply integration by parts to the integrand with
second-order derivatives,

.. math::

         -\int_\Omega (\Delta u)v dx
        = \int_\Omega\nabla u\cdot\nabla v dx - \int_{\partial\Omega}{\partial u\over
        \partial n}v ds ,

where :math:`{\partial u\over
\partial n}` is the derivative of :math:`u` in the outward normal direction at
the boundary.
The test function :math:`v` is required to vanish on the parts of the
boundary where :math:`u` is known, which in the present problem implies that
:math:`v=0` on the whole boundary :math:`\partial\Omega`.
The second term on the right-hand side of the last equation therefore
vanishes. It then follows that

.. math::

         \int_\Omega\nabla u\cdot\nabla v dx = \int_\Omega fv dx\thinspace .


This equation is supposed to hold for all :math:`v` in some function
space :math:`\hat V`. The trial function :math:`u` lies in some (possibly
different) function space :math:`V`.  We say that the last equation is
the *weak form* of the original boundary value problem consisting of
the PDE :math:`-\nabla^2u=f` and the boundary condition :math:`u=u_0`.

The proper statement of our variational problem now goes as follows:
Find :math:`u \in V` such that

.. math::

          \int_{\Omega} \nabla u \cdot \nabla v dx =
          \int_{\Omega} fv dx
          \quad \forall v \in \hat{V}.

The test and trial spaces :math:`\hat{V}` and :math:`V` are in the present
problem defined as

.. math::

            \hat{V} &= \{v \in H^1(\Omega) : v = 0 \mbox{ on } \partial\Omega\}, \\
             V      &= \{v \in H^1(\Omega) : v = u_0 \mbox{ on } \partial\Omega\}\thinspace .

In short, :math:`H^1(\Omega)` is the mathematically well-known
Sobolev space containing functions :math:`v` such that :math:`v^2` and
:math:`||\nabla v||^2` have finite integrals over :math:`\Omega`. The
solution of the underlying PDE must lie in a function space where also
the derivatives are continuous, but the Sobolev space :math:`H^1(\Omega)`
allows functions with discontinuous derivatives.  This weaker continuity
requirement of :math:`u` in the variational statement, caused by the
integration by parts, has great practical consequences when it comes to
constructing finite elements.

To solve the Poisson equation numerically, we need to transform the
continuous variational problem
to a discrete variational
problem. This is done by introducing *finite-dimensional* test and
trial spaces, often denoted as
:math:`\hat{V}_h\subset\hat{V}` and :math:`V_h\subset{V}`. The
discrete variational problem reads:
Find :math:`u_h \in V_h \subset V` such that

.. math::

          \int_{\Omega} \nabla u_h \cdot \nabla v dx =
          \int_{\Omega} fv dx
          \quad \forall v \in \hat{V}_h \subset \hat{V}\thinspace .

The choice of :math:`\hat{V}_h` and :math:`V_h` follows directly from the
kind of finite elements we want to apply in our problem. For example,
choosing the well-known linear triangular element with three nodes
implies that :math:`\hat V_h` and :math:`V_h` are the spaces of all
piecewise linear functions over a mesh of triangles, where the functions
in :math:`\hat V_h` are zero on the boundary and those in :math:`V_h`
equal :math:`u_0` on the boundary.

The mathematics literature on variational problems writes :math:`u_h` for
the solution of the discrete problem and :math:`u` for the solution of the
continuous problem. To obtain (almost) a one-to-one relationship
between the mathematical formulation of a problem and the
corresponding FEniCS program, we shall use :math:`u` for the solution of
the discrete problem and :math:`u_{e}` for the exact solution of the
continuous problem, *if* we need to explicitly distinguish
between the two.  In most cases we will introduce the PDE problem with
:math:`u` as unknown, derive a variational equation :math:`a(u,v)=L(v)` with :math:`u\in
V` and :math:`v\in \hat V`, and then simply discretize the problem by saying
that we choose finite-dimensional spaces for :math:`V` and :math:`\hat V`. This
restriction of :math:`V` implies that :math:`u` becomes a discrete finite element
function.  In practice this means that we turn our PDE problem into a
continuous variational problem, create a mesh and specify an element
type, and then let :math:`V` correspond to this mesh and element choice.
Depending upon whether :math:`V` is infinite- or finite-dimensional, :math:`u`
will be the exact or approximate solution.

It turns out to be convenient to
introduce the following unified notation for weak forms:

.. math::

        a(u, v) = L(v)\thinspace .

In the present problem we have that

.. math::

        a(u, v) &= \int_{\Omega} \nabla u \cdot \nabla v dx,
        \\
        L(v) &= \int_{\Omega} fv dx\thinspace .

From the mathematics literature, :math:`a(u,v)` is known as a *bilinear
form* and :math:`L(u)` as a *linear form*.  We shall in every problem
we solve identify the terms with the unknown :math:`u` and collect them
in :math:`a(u,v)`, and similarly collect all terms with only known
functions in :math:`L(v)`. The formulas for :math:`a` and :math:`L`
are then coded directly in the program.

To summarize, before making a FEniCS program for solving a PDE,
we must first perform two steps:

  * Turn the PDE problem into a discrete
    variational problem: Find :math:`u\in V`
    such that :math:`a(u,v) = L(v)\quad\forall v\in \hat{V}`.

  * Specify the choice of spaces (:math:`V` and :math:`\hat V`),
    i.e., the mesh and type of finite elements.

.. _tut:poisson1:impl:

Implementation (3)
------------------

The test problem so far has a general domain :math:`\Omega` and general functions
:math:`u_0` and :math:`f`. However,
we must specify :math:`\Omega`, :math:`u_0`, and :math:`f` prior to our first implementation.
It will be wise to construct a specific problem where we can easily check
that the solution is correct.
Let us choose :math:`u(x,y)=1 + x^2 + 2y^2` to be the solution of our
Poisson problem since the finite element method with linear elements
over a uniform mesh of triangular cells
should exactly reproduce a second-order polynomial
at the vertices of the cells, regardless of the size
of the elements. This property allows us to verify the code by
using very few elements and
checking that the computed and the exact solution equal to the
machine precision.
Test problems with this property will be frequently constructed throughout
the present
tutorial.
.. Should errors in the implementation arise, it is possible
.. to perform hand calculations of the intermediate steps in the finite
.. element method and compare with what the program gives.

Specifying :math:`u(x,y)=1 + x^2 + 2y^2` in the
problem from the section :ref:`tut:poisson1:varform` implies
:math:`u_0(x,y)= 1 + x^2 + 2y^2`
and :math:`f(x,y)=-6`.
We let :math:`\Omega` be the unit square for simplicity.
A FEniCS program for solving the Poisson equation in 2D
with the given choices
of :math:`u_0`, :math:`f`, and :math:`\Omega` may look as follows (the complete code can be
found in the file ``Poisson2D_D1.py``):

.. code-block:: python

        """
        FEniCS tutorial demo program:
        Poisson equation with Dirichlet conditions.
        Simplest example of computation and visualization.

        -Laplace(u) = f on the unit square.
        u = u0 on the boundary.
        u0 = u = 1 + x^2 + 2y^2, f = -6.
        """

        from dolfin import *

        # Create mesh and define function space
        mesh = UnitSquare(6, 4)
        V = FunctionSpace(mesh, 'CG', 1)

        # Define boundary conditions
        u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')

        def u0_boundary(x, on_boundary):
            return on_boundary

        bc = DirichletBC(V, u0, u0_boundary)

        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)
        f = Constant(-6.0)
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx

        # Compute solution
        problem = VariationalProblem(a, L, bc)
        u = problem.solve()

        # Plot solution and mesh
        plot(u)
        plot(mesh)

        # Dump solution to file in VTK format
        file = File('poisson.pvd')
        file << u

        # Hold plot
        interactive()

We shall now dissect this FEniCS program in detail. The program is written
in the Python programming language.  You may either take a quick look
at the `official Python tutorial <http://docs.python.org/tutorial/>`_
to pick up the basics of Python if you are unfamiliar with the language,
or you may learn enough Python as you go along with the examples in
the present tutorial. The latter strategy has proven to work for many
newcomers to FEniCS. (The requirement of using Python and an abstract
mathematical formulation of the finite element problem may seem difficult
for those who are unfamiliar with these topics.  However, the amount
of mathematics and Python that is really demanded to get you productive
with FEniCS is quited limited.  And Python is an easy-to-learn language
that you certainly will love and use far beyond FEniCS programming.)
the section :ref:`tut:appendix:pybooks` lists some relevant Python books.

The listed FEniCS program defines a finite element mesh, the discrete
function spaces :math:`V` and :math:`\hat{V}` corresponding to this
mesh and the element type, boundary conditions for :math:`u` (i.e., the
function :math:`u_0`), :math:`a(u,v)`, and :math:`L(v)`.  Thereafter,
the unknown trial function :math:`u` is computed. Then we can investigate
:math:`u` visually or analyze the computed values.

The first line in the program,

.. code-block:: python

        from dolfin import *

imports the key classes ``UnitSquare``,
``FunctionSpace``, ``Function``, and so forth, from the DOLFIN library.
All FEniCS programs for solving PDEs by the finite element method
normally start with this line. DOLFIN is a software library with efficient
and convenient C++ classes for finite element computing, and
``dolfin`` is a Python package providing access to this
C++ library from Python programs.
You can think of FEniCS as an umbrella, or project name, for a set of
computational components, where DOLFIN is one important component for
writing finite element programs. DOLFIN applies other components
in the FEniCS suite under the hood, but newcomers to FEniCS
programming do not need to care about this.

.. index:: Mesh

.. index:: DOLFIN mesh

The statement

.. code-block:: python

        mesh = UnitSquare(6, 4)

defines a uniform finite element mesh over the unit square
:math:`[0,1]\times [0,1]`. The mesh consists of *cells*, which are
triangles with straight sides. The parameters 6 and 4 tell that the
square is first divided into :math:`6\cdot 4` rectangles, and then
each rectangle is divided into two triangles. The total number of
triangles then becomes 48. The total number of vertices in this mesh is
:math:`7\cdot 5=35`.  DOLFIN offers some classes for creating meshes
over very simple geometries. For domains of more complicated shape
one needs to use a separate *preprocessor* program to create the mesh.
The FEniCS program will then read the mesh from file.

Having a mesh, we can define a discrete function space ``V`` over
this mesh:

.. index:: FunctionSpace

.. code-block:: python

        V = FunctionSpace(mesh, 'CG', 1)

The second argument reflects the type of element, while the third
argument is the degree of the basis functions on the element.

.. index:: finite element specifications

.. index:: CG finite element family

.. index:: Lagrange finite element family

Here, ``'CG'`` stands for Continuous Galerkin, implying the standard
Lagrange family of elements.  Instead of ``'CG'`` we could have written
``'Lagrange'``.  With degree 1, we simply get the standard linear
Lagrange element, which is a triangle with nodes at the three vertices.
Some finite element practitioners refer to this element as the "linear
triangle".  The computed :math:`u` will be continuous and linearly varying
in :math:`x` and :math:`y` over each cell in the mesh.  Higher-degree
polynomial approximations over each cell are trivially obtained by
increasing the third parameter in ``FunctionSpace``. Changing the second
parameter to ``'DG'`` creates a function space for discontinuous Galerkin
methods.

.. index:: TestFunction

.. index:: TrialFunction

.. index:: DirichletBC

.. index:: Dirichlet boundary conditions

In mathematics, we distinguish between the trial and test spaces :math:`V`
and :math:`\hat{V}`. The only difference in the present problem is the
boundary conditions. In FEniCS we do not specify the boundary conditions
as part of the function space, so it is sufficient to work with one
common space ``V`` for the test and trial functions in the program:

.. code-block:: python

        v = TestFunction(V)
        u = TrialFunction(V)

The next step is to specify the boundary condition: :math:`u=u_0` on
:math:`\partial\Omega`. This is done by

.. code-block:: python

        bc = DirichletBC(V, u0, u0_boundary)

where ``u0`` is an instance holding the :math:`u_0` values, and
``u0_boundary`` is a function (or object) describing whether a point
lies on the boundary where :math:`u` is specified.

Boundary conditions of the type :math:`u=u_0` are known as *Dirichlet
conditions*, and also as *essential boundary conditions* in a finite
element context.  Naturally, the name of the DOLFIN class holding the
information about Dirichlet boundary conditions is ``DirichletBC``.

.. index:: Expression

The ``u0`` variable refers to an ``Expression`` object, which is used
to represent a mathematical function. The typical construction is

.. code-block:: python

        u0 = Expression(formula)

where ``formula`` is a string containing the mathematical expression.
This formula is written with C++ syntax (the expression is automatically
turned into an efficient, compiled C++ function, see the section
:ref:`tut:app:cpp:functions` for details on the syntax). The independent
variables in the function expression are supposed to be available as
a point vector ``x``, where the first element ``x[0]`` corresponds to
the :math:`x` coordinate, the second element ``x[1]`` to the :math:`y`
coordinate, and (in a three-dimensional problem) ``x[2]`` to the :math:`z`
coordinate. With our choice of :math:`u_0(x,y)=1 + x^2 + 2y^2`, the
formula string must be written as ``1 + x[0]*x[0] + 2*x[1]*x[1]``:

.. code-block:: python

        u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')

The information about where to apply the ``u0`` function as boundary
condition is coded in a function ``boundary``:

.. index:: boundary specification (function)

.. code-block:: python

        def u0_boundary(x, on_boundary):
            return on_boundary

A function like ``u0_boundary`` for marking the boundary must return a
boolean value: ``True`` if the point ``x`` lies on the Dirichlet boundary
and ``False`` otherwise.  The argument ``on_boundary`` is ``True`` if
``x`` is on the physical boundary of the mesh, so in the present case
we can just return ``on_boundary``.  The ``u0_boundary`` function will
be called for every discrete point in the mesh, which allows us to have
boundaries where :math:`u` are known also inside the domain, if desired.

One can also omit the ``on_boundary`` argument, but in that case we need
to test on the value of the coordinates in ``x``:

.. code-block:: python

        def u0_boundary(x):
            return x[0] == 0 or x[1] == 0 or x[0] == 1 or x[1] == 1

As for the formula in ``Expression`` objects, ``x`` in the ``u0_boundary``
function represents a point in space with coordinates ``x[0]``, ``x[1]``,
etc. Comparing floating-point values using an exact match test with
``==`` is not good programming practice, because small round-off errors
in the computations of the ``x`` values could make a test ``x[0] == 1``
become false even though ``x`` lies on the boundary.  A better test is
to check for equality with a tolerance:

.. code-block:: python

        def u0_boundary(x):
            tol = 1E-15
            return abs(x[0]) < tol or \
                   abs(x[1]) < tol or \
                   abs(x[0] - 1) < tol or \
                   abs(x[1] - 1) < tol

Before defining :math:`a(u,v)` and :math:`L(v)` we have to specify the
:math:`f` function:

.. code-block:: python

        f = Expression('-6')

When :math:`f` is constant over the domain, ``f`` can be more efficiently
represented as a ``Constant`` object:

.. code-block:: python

        f = Constant(-6.0)

Now we have all the objects we need in order to specify this problem's
:math:`a(u,v)` and :math:`L(v)`:

.. code-block:: python

        a = inner(grad(u), grad(v))*dx
        L = f*v*dx

In essence, these two lines specify the PDE to be solved.  Note the very
close correspondence between the Python syntax and the mathematical
formulas :math:`\nabla u\cdot\nabla v dx` and :math:`fv dx`.  This is
a key strength of FEniCS: the formulas in the variational formulation
translate directly to very similar Python code, a feature that makes it
easy to specify PDE problems with lots of PDEs and complicated terms in
the equations.  The language used to express weak forms is called UFL
(Unified Form Language) and is an integral part of FEniCS.

Having ``a`` and ``L`` defined, and information about essential
(Dirichlet) boundary conditions in ``bc``, we can formulate a
``VariationalProblem``:

.. code-block:: python

        problem = VariationalProblem(a, L, bc)

Solving the variational problem for the solution ``u`` is just a
matter of writing

.. code-block:: python

        u = problem.solve()

Unless otherwise stated, a sparse direct solver is used to solve the
underlying linear system implied by the variational formulation. The
type of sparse direct solver depends on which linear algebra package
that is used by default. If DOLFIN is compiled with PETSc, that package
is the default linear algebra backend, otherwise it is uBLAS.  The FEniCS
distribution for Ubuntu Linux contains PETSc, and then the default solver
becomes the sparse LU solver from UMFPACK (which PETSc has an interface
to). We shall later in the section :ref:`tut:linsys` demonstrate how to
get full control of the choice of solver and any solver parameters.

The ``u`` variable refers to a finite element function, called simply a
``Function`` in FEniCS terminology.  Note that we first defined ``u``
as a ``TrialFunction`` and used it to specify ``a``.  Thereafter,
we redefined ``u`` to be a ``Function`` representing the computed
solution. This redefinition of the variable ``u`` is possible in Python
and a programming practice in FEniCS applications.

The simplest way of quickly looking at ``u`` and the mesh is to say

.. code-block:: python

        plot(u)
        plot(mesh)
        interactive()

The ``interactive()`` call is necessary for the plot to remain on the
screen. With the left, middle, and right mouse buttons you can rotate,
translate, and zoom (respectively) the plotted surface to better examine
what the solution looks like.

It is also possible to dump the computed solution to file, e.g., in the
VTK format:

.. code-block:: python

        file = File('poisson.pvd')
        file << u

The ``poisson.pvd`` file can now be loaded into any front-end to VTK, say
ParaView or VisIt. The ``plot`` function from Viper is intended for quick
examination of the solution during program development.  More in-depth
visual investigations of finite element solution will normally benefit
from using highly professional tools such as ParaView and VisIt.


.. _tut:poisson1:verify1:

Examining the Discrete Solution
-------------------------------

We know that, in the particular boundary-value problem of the section
:ref:`tut:poisson1:impl`, the computed solution :math:`u` should equal
the exact solution at the vertices of the cells.  An important extension
of our first program is therefore to examine the computed values of the
solution, which is the focus of the present section.

A finite element function like :math:`u` is expressed as a linear
combination of basis functions :math:`\phi_i`, spanning the space
:math:`V`:

.. math::

        \sum_{j=1}^N U_j \phi_j \thinspace .

By writing ``u = problem.solve()`` in the program, a linear system will
be formed from :math:`a` and :math:`L`, and this system is solved for the
:math:`U_1,\ldots,U_N` values. The :math:`U_1,\ldots,U_N` values are known

.. index:: degree of freedom

as *degrees of freedom* of :math:`u`. For Lagrange elements (and many
other element types) :math:`U_k` is simply the value of :math:`u` at
the node with global number :math:`k`.  (The nodes and cell vertices
coincide for linear Lagrange elements, while for higher-order elements
there may be additional nodes at the facets and in the interior of cells.)

Having ``u`` represented as a ``Function`` object, we can either evaluate
``u(x)`` at any vertex ``x`` in the mesh, or we can grab all the values
:math:`U_j` directly by

.. code-block:: python

        u_nodal_values = u.vector()

The result is a DOLFIN ``Vector`` object, which is basically an
encapsulation of the vector object used in the linear algebra package
that is applied to solve the linear system arising form the variational
problem.  Since we program in Python it is convenient to convert the
``Vector`` object to a standard ``numpy`` array for further processing:

.. index:: degrees of freedom array

.. index:: nodal values array

.. code-block:: python

        u_array = u_nodal_values.array()

With ``numpy`` arrays we can write "Matlab-like" code to analyze the
data. Indexing is done with square brackets: ``u_array[i]``, where the
index ``i`` always starts at ``0``.

The coordinates of the vertices in the mesh can be extracted by

.. code-block:: python

        coor = mesh.coordinates()

For a $d$-dimensional problem, ``coor`` is an :math:`M\times d` ``numpy``
array, :math:`M` being the number of vertices in the mesh. Writing out
the solution on the screen can now be done by a simple loop:

.. code-block:: python

        for i in range(len(u_array)):
            print 'u(%8g,%8g) = %g' % \
                  (coor[i][0], coor[i][1], u_array[i])

The beginning of the output looks like

.. code-block:: py

        u(       0,       0) = 1
        u(0.166667,       0) = 1.02778
        u(0.333333,       0) = 1.11111
        u(     0.5,       0) = 1.25
        u(0.666667,       0) = 1.44444
        u(0.833333,       0) = 1.69444
        u(       1,       0) = 2

For Lagrange elements of degree higher than one, the vertices and the
nodes do not coincide, and then the loop above is meaningless.

.. index:: interpolation

.. index:: interpolate

For verification purposes we want to compare the values of ``u`` at
the nodes, i.e., the values of the vector ``u_array``, with the exact
solution given by ``u0``. At each node, the difference between the
computed and exact solution should be less than a small tolerance. The
exact solution is given by the ``Expression`` object ``u0``, which we can
evaluate directly as ``u0(coor[i])`` at the vertex with global number
``i``, or as ``u0(x)`` for any spatial point.  Alternatively, we can
make a finite element field ``u_e``, representing the exact solution,
whose values at the nodes are given by the ``u0`` function. With
mathematics, :math:`u_{\mbox{e}} = \sum_{j=1}^N E_j\phi_j`, where
:math:`E_j=u_0(x_j,y_j)`, :math:`(x_j,y_j)` being the coordinates of node
number :math:`j`.  This process is known as interpolation.  FEniCS has
a function for performing the operation:

.. code-block:: python

        u_e = interpolate(u0, V)

The maximum error can now be computed as

.. code-block:: python

        u_e_array = u_e.vector().array()
        diff = abs(u_array - u_e_array)
        print 'Max error:', diff.max()

        # or more compactly:
        print 'Max error:', abs(u_e_array - u_array).max()

The value of the error should be at the level of the machine precision
(:math:`10^{-16}`).

To demonstrate the use of point evaluations of ``Function`` objects,
we write out the computed ``u`` at the center point
of the domain and compare it with the exact solution:

.. code-block:: python

        center = (0.5, 0.5)
        u_value = u(center)
        u0_value = u0(center)
        print 'numerical u at the center point:', u_value
        print 'exact     u at the center point:', u0_value

Trying a :math:`3\times 3` mesh, the output from the
previous snippet becomes

.. code-block:: py


        numerical u at the center point: [ 1.83333333]
        exact     u at the center point: [ 1.75]

The discrepancy is due to the fact that the center point is not a node
in this particular mesh, but a point in the interior of a cell,
and ``u`` varies linearly over the cell while
``u0`` is a quadratic function.

Mesh information can be gathered from the ``mesh`` object, e.g.,

  * ``mesh.num_cells()`` returns the number of cells (triangles) in the mesh,

  * ``mesh.num_vertices()`` returns the number of vertices in the mesh
    (with our choice of linear Lagrange elements this equals
    the number of nodes)

Writing ``print mesh`` dumps a short, "pretty print" description
of the mesh (``print mesh`` actually displays the result of str(mesh)`,
which defines the pretty print):

.. code-block:: py


        <Mesh of topological dimension 2 (triangles) with
        16 vertices and 18 cells, ordered>

and


.. index:: pydoc


All mesh objects are of type ``Mesh`` so typing the command
``pydoc dolfin.Mesh``
in a terminal window
will give a list of methods (i.e., functions in a class)
that can be called through any
``Mesh`` object. In fact, ``pydoc dolfin.X`` shows the
documentation of
any DOLFIN name ``X`` (at the time of this writing, some names
have missing or incomplete documentation).

We have seen how to extract the nodal values in a ``numpy`` array.
If desired, we can adjust the nodal values too. Say we want to
normalize the solution such that :math:`\max_j U_j = 1`. Then we
must divide all :math:`U_j` values
by :math:`\max_j U_j`. The following snippet performs the task:

.. code-block:: python

        max_u = u_array.max()
        u_array /= max_u
        u.vector()[:] = u_array
        print u.vector().array()

That is, we manipulate ``u_array`` as desired, and then
we insert this array into `u`'s ``Vector`` object.
The ``/=`` operator implies an
in-place modification of the object on the left-hand side: all
elements of the ``u_array`` are divided by the value ``max_u``.
Alternatively, one could write
``u_array = u_array/max_u``, which implies creating a new
array on the right-hand side and assigning this array to the
name ``u_array``.
We can equally well insert the entries of ``u_array`` into
`u`'s ``numpy`` array:

.. code-block:: python

        u.vector().array()[:] = u_array

All the code in this subsection can be found in the file ``Poisson2D_D2.py``.
.. We have commented out the ``plot`` and ``interactive`` calls in
.. this version of the program, but if you want plotting to happen, make
.. sure that ``interactive`` is called at the very end of the program.


.. _tut:poisson:membrane:

Formulating a Real Physical Problem
-----------------------------------

Perhaps you are not particularly
amazed by viewing the simple surface of :math:`u` in the
test problem from the sections :ref:`tut:poisson1:impl`
and :ref:`tut:poisson1:verify1`.
However, solving a real physical problem with a more interesting and amazing
solution on the screen
is only a matter
of specifying a more exciting domain, boundary condition, and/or
right-hand side :math:`f`.

One possible physical problem regards the deflection
:math:`D(x,y)` of an elastic circular membrane
with radius :math:`R`, subject to a localized perpendicular pressure
force, modeled as a Gaussian function.
The appropriate PDE model is

.. math::


        -T\Delta D = p(x,y)\quad\hbox{in }\Omega = \{ (x,y)\,|\, x^2+y^2\leq R\},


with

.. math::


        p(x,y) = {A\over 2\pi\sigma}\exp{\left(
        - {1\over2}\left( {x-x_0\over\sigma}\right)^2
        - {1\over2}\left( {y-y_0\over\sigma}\right)^2
        \right)}\, .


Here, :math:`T` is the tension in the membrane (constant), :math:`p` is the external
pressure load,
:math:`A` the amplitude of the pressure, :math:`(x_0,y_0)` the localization of
the Gaussian pressure function, and :math:`\sigma` the "width" of this
function. The boundary condition is :math:`D=0`.

We introduce a scaling with :math:`R` as characteristic length and
:math:`8\pi\sigma T/A` as characteristic size of :math:`D`.
(Assuming :math:`\sigma` large enough so that
:math:`p\approx\hbox{const} \sim A/(2\pi\sigma)`
in :math:`\Omega`, we can integrate an axi-symmetric version of the
equation in the radial coordinate :math:`r\in [0,R]`
and obtain :math:`D=(r^2-R^2)A/(8\pi\sigma T)`,
which for :math:`r=0` gives a rough estimate of the size of :math:`|D|`:
:math:`AR^2/(8\pi\sigma T)`.)
With this scaling we can derive the equivalent
dimensionless problem on the unit circle,

.. math::



        -\Delta w =
        4\exp{\left(
        - {1\over2}\left( {Rx-x_0\over\sigma}\right)^2
        - {1\over2}\left( {Ry-y_0\over\sigma}\right)^2
        \right)},


with :math:`w=0` on the boundary. We have that :math:`D = AR^2w/(8\pi\sigma T)`.

A mesh over the unit circle can be created
by

.. code-block:: python

        mesh = UnitCircle(n)

where ``n`` is the typical number of elements in the radial direction.
You should now be able to figure out how to modify the
``Poisson2D_D1.py`` code to solve this membrane problem.
More specifically, you are recommended to perform the following extensions:

  * initialize :math:`R`, :math:`x_0`, :math:`y_0`, :math:`\sigma`, :math:`T`, and :math:`A` in the
    beginning of the program,

  * build a string expression for :math:`p` with correct C++ syntax
    (use "printf" formatting in Python to build the expression),

  * define the ``a`` and ``L`` variables in the variational problem for
    :math:`w` and compute the solution,

  * plot the mesh, :math:`w`, and the scaled pressure function
    :math:`p` (the right-hand side of the scaled PDE),

  * write out the maximum real deflection :math:`D`
    (i.e., the maximum of the :math:`w` values times :math:`A/(8\pi\sigma T)`).

Use variable names in the program similar to the mathematical symbols
in this problem.

Choosing a small width :math:`\sigma` (say 0.01)
and a location :math:`(x_0,y_0)` toward the circular boundary
(say :math:`(0.6R\cos\theta, 0.6R\sin\theta)` for any :math:`\theta\in [0,2\pi]`),
may produce an exciting visual comparison of :math:`w` and :math:`p` that
demonstrates the very smoothed elastic response to a peak force
(or mathematically, the smoothing properties of the inverse of the
Laplace operator).
You need to experiment with the mesh resolution to get a smooth
visual representation of :math:`p`.

In the limit :math:`\sigma\rightarrow\infty`, the right-hand side function
:math:`p` approaches the constant 4,
and then the solution should be :math:`w(x,y) = 1-x^2-y^2`.
Compute the absolute value of the
difference between the exact and the numerical solution
if :math:`\sigma \geq 50` and write out the maximum difference
to provide some evidence that the implementation is correct.

You are strongly encouraged to spend some time on doing
this exercise and play around with
the plots and different mesh resolutions.
A suggested solution to the exercise
can be found in the file ``membrane1.py``.



.. code-block:: python

        """
        FEniCS program for the deflection w(x,y) of a membrane:
        -Laplace(w) = p = Gaussian function, in a unit circle,
        with w = 0 on the boundary.
        """

        from dolfin import *

        # Set pressure function:
        T = 10.0  # tension
        A = 1.0   # pressure amplitude
        R = 0.3   # radius of domain
        theta = 0.2
        x0 = 0.6*R*cos(theta)
        y0 = 0.6*R*sin(theta)
        sigma = 0.025
        #sigma = 50  # verification
        pressure = '4*exp(-0.5*(pow((%g*x[0] - %g)/%g, 2)) '\
                   '     - 0.5*(pow((%g*x[1] - %g)/%g, 2)))' % \
                   (R, x0, sigma, R, y0, sigma)

        n = 40   # approx no of elements in radial direction
        mesh = UnitCircle(n)
        V = FunctionSpace(mesh, 'CG', 1)

        # Define boundary condition w=0

        def boundary(x, on_boundary):
            return on_boundary

        bc = DirichletBC(V, Constant(0.0), boundary)

        # Define variational problem
        v = TestFunction(V)
        w = TrialFunction(V)
        p = Expression(pressure)
        a = inner(grad(w), grad(v))*dx
        L = v*p*dx

        # Compute solution
        problem = VariationalProblem(a, L, bc)
        w = problem.solve()

        # Plot solution and mesh
        plot(mesh, title='Mesh over scaled domain')
        plot(w, title='Scaled deflection')
        p = interpolate(p, V)
        plot(p, title='Scaled pressure')

        # Find maximum real deflection
        max_w = w.vector().array().max()
        max_D = A*max_w/(8*pi*sigma*T)
        print 'Maximum real deflection is', max_D

        # Verification for "flat" pressure (big sigma)
        if sigma >= 50:
            w_exact = Expression('1 - x[0]*x[0] - x[1]*x[1]')
            w_e = interpolate(w_exact, V)
            w_e_array = w_e.vector().array()
            w_array = w.vector().array()
            diff_array = abs(w_e_array - w_array)
            print 'Verification of the solution, max difference is %.4E' % \
                  diff_array.max()

            # Create finite element field over V and fill with error values
            difference = Function(V)
            difference.vector()[:] = diff_array
            #plot(difference, title='Error field for sigma=%g' % sigma)

        # Should be at the end
        interactive()





.. _tut:poisson:gradu:

Computing Derivatives
---------------------

In many Poisson and other problems the gradient of the solution is
of interest. The computation is in principle simple:
since
:math:`u = \sum_{j=1}^N U_j \phi_j`, we have that

.. math::


        \nabla u = \sum_{j=1}^N U_j \nabla \phi_j\thinspace .


Given the solution variable ``u`` in the program, ``grad(u)`` denotes
the gradient. However, the gradient of a piecewise continuous
finite element scalar field
is a discontinuous vector field
since the :math:`\phi_j` has discontinuous derivatives at the boundaries of
the cells. For example, using Lagrange elements of degree 1, :math:`u` is
linear over each cell, and the numerical :math:`\nabla u` becomes a piecewise
constant vector field. On the contrary,
the exact gradient is continuous.
For visualization and data analysis purposes
we often want the computed
gradient to be a continuous vector field. Typically,
we want each component of :math:`\nabla u` to be represented in the same
way as :math:`u` itself. To this end, we can project the components
of :math:`\nabla u` onto the
same function space as we used for :math:`u`.
This means that we solve :math:`w = \nabla u` approximately by a finite element
method, using the the same elements for the components of
:math:`w` as we used for :math:`u`. This process is known as *projection*.

.. index:: projection

Looking at the component :math:`\partial u/\partial x` of the gradient, we project
the (discrete) derivative
:math:`\sum_jU_j{\partial \phi_j/\partial x}` onto another function space
with basis :math:`\bar\phi_1,\bar\phi_2,\ldots` such that the derivative in
this space is expressed by the standard sum
:math:`\sum_j\bar U_j\bar \phi_j`, for suitable (new)
coefficients :math:`\bar U_j`.

The variational problem for :math:`w` reads: Find  :math:`w\in V^{(\mbox{g})}` such that

.. math::


        a(w, v) = L(v)\quad\forall v\in \hat{V^{(\mbox{g})}},


where

.. math::


        a(w, v) &= \int_\Omega w\cdot v dx,\\
        L(v) &= \int_\Omega \nabla u\cdot v dx\thinspace .


The function spaces :math:`V^{(\mbox{g})}` and :math:`\hat{V^{(\mbox{g})}}` (with the superscript
g denoting "gradient") are
vector versions of the function space for :math:`u`, with
boundary conditions removed (if :math:`V` is the
space we used for :math:`u`, with no restrictions
on boundary values, :math:`V^{(\mbox{g})} = \hat{V^{(\mbox{g})}} = [V]^d`, where
:math:`d` is the number of space dimensions).
For example, if we used piecewise linear functions on the mesh to
approximate :math:`u`, the variational problem for :math:`w` corresponds to
approximating each component field of :math:`w` by piecewise linear functions.

The variational problem for the vector field
:math:`w`, called ``gradu`` in the code, is easy to solve in FEniCS:

.. code-block:: python

        V_g = VectorFunctionSpace(mesh, 'CG', 1)
        v = TestFunction(V_g)
        w = TrialFunction(V_g)

        a = inner(w, v)*dx
        L = inner(grad(u), v)*dx
        problem = VariationalProblem(a, L)
        gradu = problem.solve()

        plot(gradu, title='grad(u)')

The new thing is basically that we work with a ``VectorFunctionSpace``,
since the unknown is now a vector field, instead of the
``FunctionSpace`` object for scalar fields.

The scalar component fields of the gradient
can be extracted as separated fields and, e.g., visualized:

.. code-block:: python

        gradu_x, gradu_y = gradu.split(deepcopy=True)  # extract components
        plot(gradu_x, title='x-component of grad(u)')
        plot(gradu_y, title='y-component of grad(u)')

The ``deepcopy=True`` argument signifies a *deep copy*, which is
a general term in computer science implying that a copy of the data is
returned. (The opposite, ``deepcopy=False``,
means a *shallow copy*, where
the returned objects are just pointers to the original data.)


.. index:: degrees of freedom array


.. index:: nodal values array


.. index:: degrees of freedom array (vector field)


The ``gradu_x`` and ``gradu_y`` variables behave as
``Function`` objects. In particular, we can extract the underlying
arrays of nodal values by

.. code-block:: python

        gradu_x_array = gradu_x.vector().array()
        gradu_y_array = gradu_y.vector().array()

The degrees of freedom of the ``gradu`` vector field can also be
reached by

.. code-block:: python

        gradu_array = gradu.vector().array()

but this is a flat ``numpy`` array where the degrees of freedom
for the :math:`x` component of the gradient is stored in the first part, then the
degrees of freedom of the :math:`y` component, and so on.

The program ``Poisson2D_D3.py`` extends the
code ``Poisson2D_D2.py`` from the section :ref:`tut:poisson1:verify1`
with computations and visualizations of the gradient.
Examining the arrays ``gradu_x_array``
and ``gradu_y_array``, or looking at the plots of
``gradu_x`` and
``gradu_y``, quickly reveals that
the computed ``gradu`` field does not equal the exact
gradient :math:`(2x, 4y)` in this particular test problem where :math:`u=1+x^2+2y^2`.
There are inaccuracies at the boundaries, arising from the
approximation problem for :math:`w`. Increasing the mesh resolution shows,
however, that the components of the gradient vary linearly as
:math:`2x` and :math:`4y` in
the interior of the mesh (i.e., as soon as we are one element away from
the boundary). See the section :ref:`tut:quickviz` for illustrations of
this phenomenon.


.. index:: project

.. index:: projection


Representing the gradient by the same elements as we used for the
solution is a very common step in finite element programs, so the
formation and solution of a variational problem for :math:`w` as shown above
can be replaced by a one-line call:

.. code-block:: python

        gradu = project(grad(u), VectorFunctionSpace(mesh, 'CG', 1))

The ``project`` function can take an expression involving some
finite element function in some space and project the expression onto
another space.
The applications are many, including turning discontinuous gradient
fields into continuous ones, comparing higher- and lower-order
function approximations, and transforming a higher-order finite element
solution down to a piecewise linear field, which is required by many
visualization packages.

.. _tut:poisson1:functionals:

Computing Functionals
---------------------

.. index:: functionals


After the solution :math:`u` of a PDE is computed, we often want to compute
functionals of :math:`u`, for example,

.. math::


        {1\over2}||\nabla u||^2 \equiv {1\over2}\int_\Omega \nabla u\cdot \nabla u dx,



which often reflects the some energy quantity.
Another frequently occurring functional is the error

.. math::


        ||u_{\mbox{e}}-u|| = \left(\int_\Omega (u_{\mbox{e}}-u)^2 dx\right)^{1/2},



which is of particular interest when studying convergence properties.
Sometimes the interest concerns the flux out of a part :math:`\Gamma` of
the boundary :math:`\partial\Omega`,

.. math::


        F = -\int_\Gamma p\nabla u\cdot\pmb{n} ds,



where :math:`\pmb{n}` is an outward unit normal at :math:`\Gamma` and :math:`p` is a
coefficient (see the problem in the section :ref:`tut:possion:2D:varcoeff`
for a specific example).
All these functionals are easy to compute with FEniCS, and this section
describes how it can be done.


.. index:: energy functional


*Energy Functional.* The integrand of the
energy functional
:math:`{1\over2}\int_\Omega \nabla u\cdot \nabla u dx`
is described in the UFL language in the same manner as we describe
weak forms:

.. code-block:: python

        energy = 0.5*inner(grad(u), grad(u))*dx
        E = assemble(energy, mesh=mesh)

The ``assemble`` call performs the integration.
It is possible to restrict the integration to subdomains, or parts
of the boundary, by using
a mesh function to mark the subdomains as explained in
the section :ref:`tut:poisson:mat:neumann`.
The program ``membrane2.py`` carries out the computation of
the elastic energy

.. math::


        {1\over2}||T\nabla D||^2 = {1\over2}\left({AR\over 8\pi\sigma}\right)^2
        ||\nabla w||^2


in the membrane problem from the section :ref:`tut:poisson:membrane`.


.. index:: error functional


*Convergence Estimation.* To illustrate error computations and convergence of finite element
solutions, we modify the ``Poisson2D_D3.py`` program from
the section :ref:`tut:poisson:gradu` and specify a more complicated solution,

.. math::


        u(x,y) = \sin(\omega\pi x)\sin(\omega\pi y)


on the unit square.
This choice implies :math:`f(x,y)=2\omega^2\pi^2 u(x,y)`.
With :math:`\omega` restricted to an integer
it follows that :math:`u_0=0`. We must define the
appropriate boundary conditions, the exact solution, and the :math:`f` function
in the code:

.. code-block:: python

        def boundary(x, on_boundary):
            return on_boundary

        bc = DirichletBC(V, Constant(0.0), boundary)

        omega = 1.0
        u_exact = Expression('sin(%g*pi*x[0])*sin(%g*pi*x[1])' % \
                             (omega, omega))

        f = 2*pi**2*omega**2*u_exact


The computation of
:math:`\left(\int_\Omega (u_e-u)^2 dx\right)^{1/2}`
can be done by

.. code-block:: python

        error = (u - u_exact)**2*dx
        E = sqrt(assemble(error, mesh=mesh))

However, ``u_exact`` will here be interpolated onto
the function space ``V``, i.e., the exact solution used in
the integral will vary linearly over
the cells, and not as a sine function,
if ``V`` corresponds to linear Lagrange elements.
This may yield a smaller error ``u - u_e`` than what is actually true.

More accurate representation of the exact solution is easily achieved
by interpolating the formula onto a space defined by
higher-order elements, say of third degree:

.. code-block:: python

        Ve = FunctionSpace(mesh, 'CG', degree=3)
        u_e = interpolate(u_exact, Ve)
        error = (u - u_e)**2*dx
        E = sqrt(assemble(error, mesh=mesh))


The ``u`` function will here be automatically interpolated and
represented in the
``Ve`` space. When functions in different function spaces enter
UFL expressions, they will be represented in the space of highest
order before integrations are carried out. When in doubt, we should
explicitly interpolate ``u``:

.. code-block:: python

        u_Ve = interpolate(u, Ve)
        error = (u_Ve - u_e)**2*dx


The square in the expression for ``error`` will be expanded and lead
to a lot of terms that almost cancel when the error is small, with the
potential of introducing significant round-off errors.
The function ``errornorm`` is available for avoiding this effect
by first interpolating ``u`` and ``u_exact`` to a space with
higher-order elements, then subtracting the degrees of freedom, and
then performing the integration of the error field. The usage is simple:

.. code-block:: python

        E = errornorm(u_exact, u, normtype='L2', degree=3)

At the time of this writing, ``errornorm`` does not work with
``Expression`` objects for ``u_exact``, making the function
inapplicable for most practical purposes.
Nevertheless, we can easily express the procedure explicitly:

.. code-block:: python

        def errornorm(u_exact, u, Ve):
            u_Ve = interpolate(u, Ve)
            u_e_Ve = interpolate(u_exact, Ve)
            e_Ve = Function(Ve)
            # Subtract degrees of freedom for the error field
            e_Ve.vector()[:] = u_e_Ve.vector().array() - \
                               u_Ve.vector().array()
            error = e_Ve**2*dx
            return sqrt(assemble(error, mesh=Ve.mesh()))

The ``errornorm`` procedure turns out to be identical to computing
the expression ``(u_e - u)**2*dx`` directly in
the present test case.

Sometimes it is of interest to compute the error of the
gradient field: :math:`||\nabla (u-u_{\mbox{e}})||`
(often referred to as the :math:`H^1` seminorm of the error).
Given the error field ``e_Ve`` above, we simply write

.. code-block:: python

        H1seminorm = sqrt(assemble(inner(grad(e_Ve), grad(e_Ve))*dx,
                                   mesh=mesh))


Finally, we remove all ``plot`` calls and printouts of :math:`u` values
in the original program, and
collect the computations in a function:

.. code-block:: python

        def compute(nx, ny, polynomial_degree):
            mesh = UnitSquare(nx, ny)
            V = FunctionSpace(mesh, 'CG', degree=polynomial_degree)
            ...
            Ve = FunctionSpace(mesh, 'CG', degree=3)
            E = errornorm(u_exact, u, Ve)
            return E


Calling ``compute`` for finer and finer meshes enables us to
study the convergence rate. Define the element size
:math:`h=1/n`, where :math:`n` is the number of divisions in :math:`x` and :math:`y` direction
(``nx=ny`` in the code). We perform experiments with :math:`h_0>h_1>h_2\cdots`
and compute the corresponding errors :math:`E_0, E_1, E_3` and so forth.
Assuming :math:`E_i=Ch_i^r` for unknown constants :math:`C` and :math:`r`, we can compare
two consecutive experiments, :math:`E_i=Ch_i^r` and :math:`E_{i-1}=Ch_{i-1}^r`,
and solve for :math:`r`:

.. math::


        r = {\ln(E_i/E_{i-1})\over\ln (h_i/h_{i-1})}\thinspace .


The :math:`r` values should approach the expected convergence
rate ``degree+1`` as :math:`i` increases.

The procedure above can easily be turned into Python code:

.. code-block:: python

        import sys
        degree = int(sys.argv[1])  # read degree as 1st command-line arg
        h = []  # element sizes
        E = []  # errors
        for nx in [4, 8, 16, 32, 64, 128, 264]:
            h.append(1.0/nx)
            E.append(compute(nx, nx, degree))

        # Convergence rates
        from math import log as ln  # (log is a dolfin name too - and logg :-)
        for i in range(1, len(E)):
            r = ln(E[i]/E[i-1])/ln(h[i]/h[i-1])
            print 'h=%10.2E r=.2f'  (h[i], r)

The resulting program has the name ``Poisson2D_D4.py``
and computes error norms in various ways. Running this
program for elements of first degree and :math:`\omega=1` yields the output

.. code-block:: py


        h=1.25E-01 E=3.25E-02 r=1.83
        h=6.25E-02 E=8.37E-03 r=1.96
        h=3.12E-02 E=2.11E-03 r=1.99
        h=1.56E-02 E=5.29E-04 r=2.00
        h=7.81E-03 E=1.32E-04 r=2.00
        h=3.79E-03 E=3.11E-05 r=2.00

That is, we approach the expected second-order convergence of linear
Lagrange elements as the meshes become sufficiently fine.

Running the program for second-degree elements results in the expected
value :math:`r=3`,

.. code-block:: py


        h=1.25E-01 E=5.66E-04 r=3.09
        h=6.25E-02 E=6.93E-05 r=3.03
        h=3.12E-02 E=8.62E-06 r=3.01
        h=1.56E-02 E=1.08E-06 r=3.00
        h=7.81E-03 E=1.34E-07 r=3.00
        h=3.79E-03 E=1.53E-08 r=3.00

However, using ``(u - u_exact)**2`` for the error computation, which
implies interpolating ``u_exact`` onto the same space as ``u``,
results in :math:`r=4` (!). This is an example where it is important to
interpolate ``u_exact`` to a higher-order space (polynomials of
degree 3 are sufficient here) to avoid computing a too optimistic
convergence rate. Looking at the error in the degrees of
freedom (``u.vector().array()``) reveals a convergence rate of :math:`r=4`
for second-degree elements. For elements of polynomial degree 3
all the rates are
:math:`r=4`, regardless of whether we choose a "fine" space
``Ve`` with polynomials of degree 3 or 5.


Running the program for third-degree elements results in the
expected value :math:`r=4`:

.. code-block:: py


        h=  1.25E-01 r=4.09
        h=  6.25E-02 r=4.03
        h=  3.12E-02 r=4.01
        h=  1.56E-02 r=4.00
        h=  7.81E-03 r=4.00

Checking convergence rates is the next best method for verifying PDE codes
(the best being exact recovery of a solution as in the section :ref:`tut:poisson1:verify1` and many other places in this tutorial).


.. index:: flux functional


*Flux Functionals.* To compute flux integrals like
\int_\Gamma p\nabla u\cdot\pmb{n} ds
we need to define the :math:`\pmb{n}` vector, referred to as *facet normal*
in FEniCS. If :math:`\Gamma` is the complete boundary we can perform
the flux computation by

.. code-block:: python

        n = FacetNormal(mesh)
        flux = -p*inner(grad(u), n)*ds
        total_flux = assemble(flux)

It is possible to restrict the integration to a part of the boundary
using a mesh function to mark the relevant part, as
explained in the section :ref:`tut:poisson:mat:neumann`. Assuming that the
part corresponds to subdomain number ``n``, the relevant form for the
flux is ``-p*inner(grad(u), n)*ds(n)``.


.. _tut:quickviz:

Quick Visualization with VTK
----------------------------

.. index:: visualization

.. index:: Viper

.. index:: VTK


As we go along with examples it is fun to play around with
``plot`` commands and visualize what is computed. This section explains
some useful visualization features.

The ``plot(u)`` command launches a FEniCS component called Viper, which
applies the VTK package to visualize finite element functions.
Viper is not a full-fledged, easy-to-use front-end to VTK (like ParaView
or VisIt), but rather a thin layer on top of VTK's Python interface,
allowing us to quickly visualize a DOLFIN function or mesh, or data in
plain Numerical Python arrays, within a Python program.
Viper is ideal for debugging, teaching, and initial scientific investigations.
The visualization can be interactive, or you can steer and automate it
through program statements.
More advanced and professional visualizations are usually better done with
advanced tools like Mayavi2, ParaView, or VisIt.

We have made a program ``membrane1v.py`` for the membrane deflection
problem in the section :ref:`tut:poisson:membrane` and added various
demonstrations of Viper capabilities. You are encouraged to play around with
``membrane1v.py`` and modify the code as you read about various features.
The ``membrane1v.py`` program solves the two-dimensional Poisson
equation for a scalar field ``w`` (the membrane deflection).


.. index:: plot


The ``plot`` function can take additional arguments, such as
a title of the plot, or a specification of a wireframe plot (elevated mesh)
instead of a colored surface plot:

.. code-block:: python

        plot(mesh, title='Finite element mesh')
        plot(w, wireframe=True, title='solution')


The three mouse buttons can be used to rotate, translate, and zoom
the surface.
Pressing ``h`` in the plot window makes a printout of several
key bindings that are available in such windows. For example,
pressing ``m`` in the mesh plot window
dumps the plot of the mesh to an Encapsulated PostScript (``.eps``)
file, while pressing ``i`` saves the plot in PNG format.
All plotfile names are automatically generated as ``simulationX.eps``,
where ``X`` is a counter ``0000``, ``0001``, ``0002``, etc.,
being increased every time a new plot file in that format
is generated (the extension
of PNG files is ``.png`` instead of ``.eps``).
Pressing ``'o'`` adds a red outline of a bounding box around the domain.

One can alternatively control the visualization from the program code
directly. This is done through a ``Viper`` object returned from
the ``plot`` command. Let us grab this object and use it to
1) tilt the camera :math:`-65` degrees in latitude direction, 2) add
:math:`x` and :math:`y` axes, 3) change the default name of the plot files (generated
by typing ``m`` and ``i`` in the plot window),
4) change the color scale, and 5) write the plot
to a PNG and an EPS file. Here is the code:

.. code-block:: python

        viz_w = plot(w,
                    wireframe=False,
                    title='Scaled membrane deflection',
                    rescale=False,
                    axes=True,              # include axes
                    basename='deflection',  # default plotfile name
                    )

        viz_w.elevate(-65) # tilt camera -65 degrees (latitude dir)
        viz_w.set_min_max(0, 0.5*max_w)  # color scale
        viz_w.update(w)    # bring settings above into action
        viz_w.write_png('deflection.png')
        viz_w.write_ps('deflection', format='eps')

The ``format`` argument in the latter line can also take the values
``'ps'`` for a standard PostScript file and ``'pdf'`` for
a PDF file.
Note the necessity of the ``viz_w.update(w)`` call -- without it we will
not see the effects of tilting the camera and changing the color scale.
Figure :ref:`tut:poisson:2D:fig1` shows the resulting scalar surface.

.. parameters['plot_filename_prefix'] = 'hello' # does not work



.. _tut:poisson:2D:fig1:

.. figure:: eps/membrane_waxis.png
   :width: 400

   Plot of the deflection of a membrane



.. _tut:poisson1:DN:

Combining Dirichlet and Neumann Conditions
------------------------------------------

Let us make a slight extension of our two-dimensional Poisson problem
from the section :ref:`tut:poisson1:bvp`
and add a Neumann boundary condition. The domain is still
the unit square, but now we set the Dirichlet condition
:math:`u=u_0` at the left and right sides,
:math:`x=0` and :math:`x=1`, while the Neumann condition

.. math::


        -{\partial u\over\partial n}=g


is applied to the remaining
sides :math:`y=0` and :math:`y=1`.
The Neumann condition is also known as a *natural boundary condition*
(in contrast to an essential boundary condition).

.. index:: Neumann boundary conditions


Let :math:`\Gamma_D` and :math:`\Gamma_N`
denote the parts of :math:`\partial\Omega` where the Dirichlet and Neumann
conditions apply, respectively.
The complete boundary-value problem can be written as

.. math::


            - \Delta u =& f \mbox{ in } \Omega,  \\
            u =& u_0 \mbox{ on } \Gamma_D,       \\
            - {\partial u\over\partial n}  &=  g \mbox{ on } \Gamma_N  \thinspace .


Again we choose :math:`u=1+x^2 + 2y^2` as the exact solution and adjust :math:`f`, :math:`g`, and
:math:`u_0` accordingly:

.. math::


        f &= -6,\\
        g &= \left\lbrace\begin{array}{ll}
        -4, & y=1\\
        0,  & y=0
        \end{array}\right.\\
        u_0 =& 1 + x^2 + 2y^2\thinspace .


For ease of programming we may introduce a :math:`g` function defined over the whole
of :math:`\Omega` such that :math:`g` takes on the right values at :math:`y=0` and
:math:`y=1`. One possible extension is

.. math::


        g(x,y) = -4y\thinspace .



The first task is to derive the variational problem. This time we cannot
omit the boundary term arising from the integration by parts, because
:math:`v` is only zero at the :math:`\Gamma_D`. We have

.. math::


         -\int_\Omega (\Delta u)v dx
        = \int_\Omega\nabla u\cdot\nabla v dx - \int_{\partial\Omega}{\partial u\over
        \partial n}v ds,


and since :math:`v=0` on :math:`\Gamma_D`,

.. math::


        - \int_{\partial\Omega}{\partial u\over
        \partial n}v ds
        =
        - \int_{\Gamma_N}{\partial u\over
        \partial n}v ds
        = \int_{\Gamma_N}gv ds,


by applying the boundary condition at :math:`\Gamma_N`.
The resulting weak form reads

.. math::


        \int_{\Omega} \nabla u \cdot \nabla v dx +
        \int_{\Gamma_N} gv ds
        = \int_{\Omega} fv dx\thinspace .



Expressing this equation
in the standard notation :math:`a(u,v)=L(v)` is straightforward with

.. math::


        a(u, v) &= \int_{\Omega} \nabla u \cdot \nabla v dx,
        \\
        L(v) &= \int_{\Omega} fv dx -
        \int_{\Gamma_N} gv ds\thinspace .



How does the Neumann condition impact the implementation?
The code in the file ``Poisson2D_D2.py`` remains almost the same.
Only two adjustments are necessary:

  * The function describing the boundary where Dirichlet conditions
    apply must be modified.

  * The new boundary term must be added to the expression in ``L``.

Step 1 can be coded as

.. code-block:: python

        def Dirichlet_boundary(x, on_boundary):
            if on_boundary:
                if x[0] == 0 or x[0] == 1:
                    return True
                else:
                    return False
            else:
                return False

A more compact implementation reads

.. code-block:: python

        def Dirichlet_boundary(x, on_boundary):
            return on_boundary and (x[0] == 0 or x[0] == 1)

As pointed out already in the section :ref:`tut:poisson1:impl`,
testing for an exact match of real numbers is
not good programming practice so we introduce a tolerance in the test:

.. code-block:: python

        def Dirichlet_boundary(x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and \
                   (abs(x[0]) < tol or abs(x[0] - 1) < tol)

We may also split the boundary functions into two separate pieces, one
for each part of the boundary:

.. code-block:: python

        tol = 1E-14
        def Dirichlet_boundary0(x, on_boundary):
            return on_boundary and abs(x[0]) < tol

        def Dirichlet_boundary1(x, on_boundary):
            return on_boundary and abs(x[0] - 1) < tol

        bc0 = DirichletBC(V, Constant(0), Dirichlet_boundary0)
        bc1 = DirichletBC(V, Constant(1), Dirichlet_boundary1)
        bc = [bc0, bc1]






The second adjustment of our program concerns the definition of ``L``,
where we have to add a boundary integral and a definition of the :math:`g`
function to be integrated:

.. code-block:: python

        g = Expression('-4*x[1]')
        L = f*v*dx - g*v*ds

The ``ds`` variable implies a boundary integral, while ``dx``
implies an integral over the domain :math:`\Omega`.
No more modifications are necessary. Running the resulting program,
found in the file ``Poisson2D_DN1.py``, shows a
successful verification --
:math:`u` equals the exact solution at all the nodes, regardless of
how many elements we use.

.. _tut:poisson:multiple:Dirichlet:

Multiple Dirichlet Conditions
-----------------------------

The PDE problem from the previous section applies a function :math:`u_0(x,y)`
for setting Dirichlet conditions at two parts of the boundary.
Having a single function to set multiple Dirichlet conditions is
seldom possible. The more general case is to have :math:`m` functions for
setting Dirichlet conditions at :math:`m` parts of the boundary.
The purpose of this section is to explain how such multiple conditions
are treated in FEniCS programs.

Let us
return to the case from the section :ref:`tut:poisson1:DN`
and define two separate functions for
the two Dirichlet conditions:

.. math::


            - \Delta u &= -6 \mbox{ in } \Omega, \\
            u &= u_L \mbox{ on } \Gamma_0, \\
            u &= u_R \mbox{ on } \Gamma_1, \\
            - {\partial u\over\partial n}  &=  g \mbox{ on } \Gamma_N \thinspace .


Here, :math:`\Gamma_0` is the boundary :math:`x=0`, while
:math:`\Gamma_1` corresponds to the boundary :math:`x=1`.
We have that :math:`u_L = 1 + 2y^2`, :math:`u_R = 2 + 2y^2`, and :math:`g=-4y`.
For the left boundary :math:`\Gamma_0` we
define
the usual triple of a function for the boundary value,
a function for defining
the boundary of interest, and a ``DirichletBC`` object:

.. code-block:: python

        u_L = Expression('1 + 2*x[1]*x[1]')

        def left_boundary(x, on_nboundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[0]) < tol

        Gamma_0 = DirichletBC(V, u_L, left_boundary)

For the boundary :math:`x=1` we define a similar code:

.. code-block:: python

        u_R = Expression('2 + 2*x[1]*x[1]')

        def right_boundary(x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[0] - 1) < tol

        Gamma_1 = DirichletBC(V, u_R, right_boundary)

The various essential conditions are then collected in a list
and passed onto our problem object of type ``VariationalProblem``:

.. code-block:: python

        bc = [Gamma_0, Gamma_1]
        ...
        problem = VariationalProblem(a, L, bc)


If the :math:`u` values are constant at a part of the boundary, we may use
a simple ``Constant`` object instead of an ``Expression`` object.

The file ``Poisson2D_DN2.py`` contains a complete program which
demonstrates the constructions above.
An extended example with multiple Neumann conditions would have
been quite natural now, but this requires marking various parts
of the boundary using the mesh function concept and is therefore
left to the section :ref:`tut:poisson:mat:neumann`.


.. _tut:poisson1:linalg:

A Linear Algebra Formulation
----------------------------

Given :math:`a(u,v)=L(v)`, the discrete solution :math:`u` is computed by
inserting :math:`u=\sum_{j=1}^N U_j \phi_j` into :math:`a(u,v)` and demanding
:math:`a(u,v)=L(v)` to be fulfilled for :math:`N` test functions
:math:`\hat\phi_1,\ldots,\hat\phi_N`. This implies

.. math::


        \sum_{j=1}^N a(\phi_j,\hat\phi_i) U_j = L(\hat\phi_i),\quad i=1,\ldots,N,


which is nothing but a linear system,

.. math::


          AU = b,


where the entries in :math:`A` and :math:`b` are given by

.. math::


          A_{ij} &= a(\phi_j, \hat{\phi}_i), \\
          b_i &= L(\hat\phi_i)\thinspace .




.. index:: assemble


.. index:: linear systems (in FEniCS)


.. index:: assembly of linear systems


The examples so far have constructed a ``VariationalProblem`` object
and called its ``solve`` method to compute the solution
``u``.
The ``VariationalProblem`` object creates a linear system
:math:`AU=b` and calls an appropriate solution method for such systems.
An alternative is dropping the use of a ``VariationalProblem``
object and instead asking
FEniCS to create the matrix :math:`A`
and right-hand side :math:`b`, and then solve for the
solution vector :math:`U` of the linear system.
The relevant statements read

.. code-block:: python

        A = assemble(a)
        b = assemble(L)
        bc.apply(A, b)
        u = Function(V)
        solve(A, u.vector(), b)

The variables ``a`` and ``L`` are as before, i.e., ``a`` refers to the
bilinear form involving a ``TrialFunction`` object (say ``u``)
and a ``TestFunction`` object (``v``), and ``L`` involves a
``TestFunction`` object (``v``). From ``a`` and ``L``,
the ``assemble`` function can
compute the matrix elements :math:`A_{i,j}` and the vector elements :math:`b_i`.

The matrix :math:`A` and vector :math:`b` are first assembled without incorporating
essential (Dirichlet) boundary conditions. Thereafter, the
``bc.apply(A, b)`` call performs the necessary modifications to
the linear system. The first three statements above can alternatively
be carried out by

.. code-block:: python

        A, b = assemble_system(a, L, bc)

The essential boundary conditions are
now applied to the element matrices and vectors prior to assembly.

.. index:: assemble_system


When we have multiple Dirichlet conditions stored in a list ``bc``,
as explained in
the section :ref:`tut:poisson:multiple:Dirichlet`, we must apply
each condition in ``bc`` to the system:

.. code-block:: python

        # bc is a list of DirichletBC objects
        for condition in bc:
            condition.apply(A, b)

Alternatively, we can make the call

.. code-block:: python

        A, b = assemble_system(a, L, bc)


Note that the solution ``u`` is, as before, a ``Function`` object.
The degrees of freedom, :math:`U=A^{-1}b`, are filled
into `u`'s ``Vector`` object (``u.vector()``)
by the ``solve`` function.

The object ``A`` is of type ``Matrix``, while ``b`` and
``u.vector()`` are of type ``Vector``. We may convert the
matrix and vector data to ``numpy`` arrays by calling the
``array()`` method as shown before. If you wonder how essential
boundary conditions are incorporated in the linear system, you can
print out ``A`` and ``b`` before and after the
``bc.apply(A, b)`` call:

.. code-block:: python

        if mesh.num_cells() < 16:  # print for small meshes only
            print A.array()
            print b.array()
        bc.apply(A, b)
        if mesh.num_cells() < 16:
            print A.array()
            print b.array()

You will see that ``A`` is modified in a symmetric way:
for each degree of freedom that is known, the corresponding row
and column is zero'ed out and 1 is placed on the main diagonal.
The right-hand side ``b`` is modified accordingly (the column times
the value of the degree of freedom is subtracted from ``b``, and
then the corresponding entry in ``b`` is replaced by the known value
of the degree of freedom).


.. index:: File


Sometimes it can be handy to transfer the linear system to Matlab or Octave
for further analysis, e.g., computation of eigenvalues of :math:`A`.
This is easily done by opening
a ``File`` object with a filename extension ``.m`` and dump
the ``Matrix`` and ``Vector`` objects as follows:

.. code-block:: python

        mfile = File('A.m'); mfile << A
        mfile = File('b.m'); mfile << b

The data files ``A.m`` and ``b.m`` can be loaded directly into
Matlab or Octave.

The complete code where our Poisson problem is solved by forming
the linear system :math:`AU=b` explicitly, is stored in the files
``Poisson2D_DN_la1.py`` (one common Dirichlet condition) and
``Poisson2D_DN_la2.py`` (two separate Dirichlet conditions).

Creating the linear system
explicitly in the user's program, as an alternative to
using a ``VariationalProblem`` object, can have some advantages in more
advanced problem settings. For example, :math:`A` may be constant throughout
a time-dependent simulation, so we can avoid recalculating :math:`A` at
every time level and save a significant amount of simulation time. The sections :ref:`tut:timedep:diffusion1:impl` and
:ref:`tut:timedep:diffusion1:noassemble` deal with this topic in detail.

.. In other problems, we may divide the variational
.. problem and linear system into different terms, say :math:`A=M + {\Delta t} K`,
.. where :math:`M` is a matrix arising from a term like :math:`\partial u/\partial t`,
.. :math:`K` is a term corresponding to a Laplace operator, and :math:`{\Delta t}` is
.. a time discretization parameter. When :math:`{\Delta t}` is changed in time,
.. we can efficiently recompute :math:`A = M + {\Delta t} K` without
.. reassembling the constant matrices :math:`M` and :math:`K`. This strategy may
.. speed up simulations significantly.


.. _tut:possion:2D:varcoeff:

A Variable-Coefficient Poisson Problem
--------------------------------------

.. index:: Poisson's equation with variable coefficient


Suppose we have a variable coefficient :math:`p(x,y)` in the Laplace operator,
as in the boundary-value problem

.. math::


          \begin{array}{rcll}
            - \nabla\cdot \left\lbrack
        p(x,y)\nabla u(x,y)\right\rbrack  &=  f(x,y) &\mbox{in } \Omega, \\
            u(x,y)  &=  u_0(x,y) &\mbox{on } \partial\Omega\thinspace .
          \end{array}


We shall quickly demonstrate that this simple extension of our model
problem only requires an equally simple extension of the FEniCS program.

Let us continue to use our favorite solution :math:`u(x,y)=1+x^2+2y^2` and
then prescribe :math:`p(x,y)=x+y`. It follows that
:math:`u_0(x,y) = 1 + x^2 + 2y^2` and :math:`f(x,y)=-8x-10y`.

What are the modifications we need to do in the ``Poisson2D_D2.py`` program
from the section :ref:`tut:poisson1:verify1`?

  * ``f`` must be an ``Expression`` since it is no longer a constant,

  * a new ``Expression`` `p` must be defined for the variable coefficient,

  * the variational problem is slightly changed.

First we address the modified variational problem. Multiplying
the PDE by a test function :math:`v` and
integrating by parts now results
in

.. math::


        \int_\Omega p\nabla u\cdot\nabla v dx -
        \int_{\partial\Omega} p{\partial u\over
        \partial n}v ds = \int_\Omega fv dx\thinspace .


The function spaces for :math:`u` and :math:`v` are the same as in
the section :ref:`tut:poisson1:varform`, implying that the boundary integral
vanishes since :math:`v=0` on :math:`\partial\Omega` where we have Dirichlet conditions.
The weak form :math:`a(u,v)=L(v)` then has

.. math::


        a(u,v) &= \int_\Omega p\nabla u\cdot\nabla v dx,\\
        L(v) &= \int_\Omega fv dx\thinspace .


In the code from the section :ref:`tut:poisson1:impl` we must replace

.. code-block:: python

        a = inner(grad(u), grad(v))*dx

by

.. code-block:: python

        a = p*inner(grad(u), grad(v))*dx

The definitions of ``p`` and ``f`` read

.. code-block:: python

        p = Expression('x[0] + x[1]')
        f = Expression('-8*x[0] - 10*x[1]')

No additional modifications are necessary. The complete code can be
found in in the file ``Poisson2D_Dvc.py``. You can run it and confirm
that it recovers the exact :math:`u` at the nodes.

The flux :math:`-p\nabla u` may be of particular interest in variable-coefficient
Poisson
problems. As explained in the section :ref:`tut:poisson:gradu`,
we normally want the piecewise discontinuous flux or gradient
to be approximated by a continuous vector field, using the same elements
as used for the numerical solution :math:`u`. The approximation now consists of
solving :math:`w = -p\nabla u` by a finite element method:
find :math:`w\in V^{(\mbox{g})}` such that

.. math::


        a(w, v) = L(v)\quad\forall v\in \hat{V^{(\mbox{g})}},


where

.. math::


        a(w, v) &= \int_\Omega w\cdot v dx,\\
        L(v) &= \int_\Omega (-p \nabla u)\cdot v dx\thinspace .


This problem is identical to the one in the section :ref:`tut:poisson:gradu`,
except that :math:`p` enters the integral in :math:`L`.

The relevant Python statements for computing the flux field take the form

.. code-block:: python

        V_g = VectorFunctionSpace(mesh, 'CG', 1)
        v = TestFunction(V_g)
        w = TrialFunction(V_g)

        a = inner(w, v)*dx
        L = inner(-p*grad(u), v)*dx
        problem = VariationalProblem(a, L)
        flux = problem.solve()

The convenience function ``project`` was made to condense the frequently
occurring statements above:

.. code-block:: python

        flux = project(-p*grad(u),
                       VectorFunctionSpace(mesh, 'CG', 1))


Plotting the flux vector field is naturally as easy as plotting
the gradient (see the section :ref:`tut:poisson:gradu`):

.. code-block:: python

        plot(flux, title='flux field')

        flux_x, flux_y = flux.split(deepcopy=True)  # extract components
        plot(flux_x, title='x-component of flux (-p*grad(u))')
        plot(flux_y, title='y-component of flux (-p*grad(u))')


Data analysis of the nodal values of the flux field may conveniently
apply the underlying ``numpy`` arrays:

.. code-block:: python

        flux_x_array = flux_x.vector().array()
        flux_y_array = flux_y.vector().array()


The program ``Poisson2D_Dvc.py`` contains in addition some plots,
including a curve plot
comparing ``flux_x`` and the exact counterpart along the line :math:`y=1/2`.
The associated programming details related to this visualization
are explained in the section :ref:`tut:structviz`.

.. _tut:structviz:

Visualization of Structured Mesh Data
-------------------------------------

.. index:: structured mesh


.. index:: visualization, structured mesh


When finite element computations are done on a structured rectangular
mesh, maybe with uniform partitioning, VTK-based tools for completely
unstructured 2D/3D meshes are not required.  Instead we can use
many alternative high-quality
visualization tools for structured data, like the data appearing in
finite difference simulations and image analysis.  We shall
demonstrate the potential of such tools and how they allow for
more tailored and flexible visualization and data analysis.

A necessary first step is to transform our ``mesh`` object to an
object representing a rectangle with equally-shaped *rectangular*
cells.  The Python package ``scitools`` has this type of structure,
called a ``UniformBoxGrid``. The second step is to transform the
one-dimensional array of nodal values to a two-dimensional array
holding the values at the corners of the cells in the structured
grid. In such grids, we want to access a value by its :math:`i` and :math:`j`
indices, :math:`i` counting cells in the :math:`x` direction, and :math:`j` counting
cells in the :math:`y` direction.  This transformation is in principle
straightforward, yet it frequently leads to obscure indexing
errors. The ``BoxField`` object in ``scitools`` takes conveniently
care of the details of the transformation.  With a ``BoxField``
defined on a ``UniformBoxGrid`` it is very easy to call up more
standard plotting packages to visualize the solution along lines in
the domain or as 2D contours or lifted surfaces.

Let us go back to the ``Poisson2D_Dvc.py`` code from
the section :ref:`tut:possion:2D:varcoeff` and map ``u`` onto a
``BoxField`` object:

.. code-block:: python

        from scitools.BoxField import *
        u2 = u if u.ufl_element().degree() == 1 else \
             interpolate(u, FunctionSpace(mesh, 'CG', 1))
        u_box = dolfin_function2BoxField(u2, mesh, (nx,ny), uniform_mesh=True)

Note that the function ``dolfin_function2BoxField`` can only work with
finite element fields with *linear* (degree 1) elements, so for
higher-degree elements we here simply interpolate the solution onto
a mesh with linear elements. We could also project ``u`` or
interpolate/project onto a finer mesh in the higher-degree case.
Such transformations to linear finite element fields
are very often needed when calling up plotting packages or data analysis tools.
The ``u.ufl_element()`` method returns an object holding the element
type, and this object has a method ``degree()`` for returning the
element degree as an integer.
The parameters ``nx`` and ``ny`` are the number of divisions in each space
direction that were used when calling ``UnitSquare`` to make the
``mesh`` object.
The result ``u_box`` is a ``BoxField``
object that supports "finite difference" indexing and an underlying
grid suitable for ``numpy`` operations on 2D data.
Also 1D and 3D functions (with linear elements) in DOLFIN can be turned
into ``BoxField`` objects for plotting and analysis.

The ability to access a finite element field in the way one can access
a finite difference-type of field is handy in many occasions, including
visualization and data analysis.
Here is an example of writing out the coordinates and the field value
at a grid point with indices ``i`` and ``j`` (going from 0 to
``nx`` and ``ny``, respectively, from lower left to upper right corner):

.. code-block:: python

        i = nx; j = ny   # upper right corner
        print 'u(%g,%g)=%g' % (u_box.grid.coor[X][i],
                               u_box.grid.coor[Y][j],
                               u_box.values[i,j])

For instance,
the :math:`x` coordinates are reached by ``u_box.grid.coor[X]``, where
``X`` is an integer (0) imported from ``scitools.BoxField``.
The ``grid`` attribute is an instance of class ``UniformBoxGrid``.

Many plotting programs can be used to visualize the data in
``u_box``.  Matplotlib is now a very popular plotting program in
the Python world and could be used to make contour plots of
``u_box``. However, other programs like Gnuplot, VTK, and Matlab have better
support for surface plots. Our choice in this tutorial is to use the
Python package ``scitools.easyviz``, which offers a uniform
Matlab-like syntax as interface to various plotting packages such as Gnuplot,
Matplotlib, VTK, OpenDX, Matlab, and others. With ``scitools.easyviz`` we
write one set of statements, close to what one would do in Matlab or
Octave, and then it is easy to switch between different plotting
programs, at a later stage, through a command-line option, a line in a
configuration file, or an import statement in the program.  By
default, ``scitools.easyviz`` employs Gnuplot as plotting program,
and this is a highly relevant choice for scalar fields over two-dimensional,
structured meshes, or for curve plots along lines through the domain.


.. index:: contour plot


A contour plot is made by the following ``scitools.easyviz`` command:

.. code-block:: python

        from scitools.easyviz import contour, title, hardcopy
        contour(u_box.grid.coorv[X], u_box.grid.coorv[Y], u_box.values,
                5, clabels='on')
        title('Contour plot of u')
        hardcopy('u_contours.eps')

        # or more compact syntax:
        contour(u_box.grid.coorv[X], u_box.grid.coorv[Y], u_box.values,
                5, clabels='on',
                hardcopy='u_contours.eps', title='Contour plot of u')

The resulting plot can be viewed in Figure :ref:`tut:poisson:2D:fig2`.
The ``contour`` function needs arrays with the :math:`x` and :math:`y`
coordinates expanded to 2D arrays (in the same way as demanded when
making vectorized ``numpy`` calculations of arithmetic expressions
over all grid points).  The correctly expanded arrays are stored in
``grid.coorv``.  The above call to ``contour`` creates 5 equally
spaced contour lines, and with ``clabels='on'`` the contour values can
be seen in the plot.

Other functions for visualizing 2D scalar fields are ``surf`` and
``mesh`` as known from Matlab. Because the ``from dolfin import *``
statement imports several names that are also present
in ``scitools.easyviz`` (e.g., ``plot``, ``mesh``, and
``figure``), we use functions from the latter package through a
module prefix ``ev`` (for \underline{e}asy\underline{v}iz) from now on:

.. code-block:: python

        import scitools.easyviz as ev
        ev.figure()
        ev.surf(u_box.grid.coorv[X], u_box.grid.coorv[Y], u_box.values,
                shading='interp', colorbar='on',
                title='surf plot of u', hardcopy='u_surf.eps')

        ev.figure()
        ev.mesh(u_box.grid.coorv[X], u_box.grid.coorv[Y], u_box.values,
                title='mesh plot of u', hardcopy='u_mesh.eps')

Figure :ref:`tut:poisson:2D:fig3` exemplifies the surfaces arising from
the two plotting commands above.
You can type
``pydoc scitools.easyviz`` in a terminal window
to get a full tutorial.

A handy feature of ``BoxField`` is the ability to give a start point
in the grid and a direction, and then extract the field and corresponding
coordinates along the nearest grid
line. In 3D fields
one can also extract data in a plane.
Say we
want to plot :math:`u` along the line :math:`y=1/2` in the grid. The grid points,
``x``, and the
:math:`u` values along this line, ``uval``, are extracted by

.. code-block:: python

        start = (0, 0.5)
        x, uval, y_fixed, snapped = u_box.gridline(start, direction=X)

The variable ``snapped`` is true if the line had to be snapped onto a
gridline and in that case ``y_fixed`` holds the snapped
(altered) :math:`y` value.
Plotting :math:`u` versus the :math:`x` coordinate along this line, using
``scitools.easyviz``, is now a matter of

.. code-block:: python

        ev.figure()  # new plot window
        ev.plot(x, uval, 'r-')  # 'r--: red solid line
        ev.title('Solution')
        ev.legend('finite element solution')

        # or more compactly:
        ev.plot(x, uval, 'r-', title='Solution',
                legend='finite element solution')


A more exciting plot compares the projected numerical flux in
:math:`x` direction along the
line :math:`y=1/2` with the exact flux:

.. code-block:: python

        ev.figure()
        flux2_x = flux_x if flux_x.ufl_element().degree() == 1 else \
            interpolate(flux_x, FunctionSpace(mesh, 'CG', 1))
        flux_x_box = dolfin_function2BoxField(flux2_x, mesh, (nx,ny),
                                              uniform_mesh=True)
        x, fluxval, y_fixed, snapped = \
              flux_x_box.gridline(start, direction=X)
        y = y_fixed
        flux_x_exact = -(x + y)*2*x
        ev.plot(x, fluxval, 'r-',
                x, flux_x_exact, 'b-',
                legend=('numerical (projected) flux', 'exact flux'),
                title='Flux in x-direction (at y=%g)' % y_fixed,
                hardcopy='flux.eps')

As seen from Figure :ref:`tut:poisson:2D:fig2`, the numerical flux
is accurate except in the elements closest to the boundaries.




.. figure:: eps/Poisson2D_Dvc_contour1.png
   :width: 400

   Examples on plots created by transforming the finite element field to a field on a uniform, structured 2D grid: contour plot of the solution



.. _tut:poisson:2D:fig2:

.. figure:: eps/Poisson2D_Dvc_flux_x.png
   :width: 400

   Examples on plots created by transforming the finite element field to a field on a uniform, structured 2D grid: curve plot of the exact flux :math:`-p\partial u/\partial x` against the corresponding projected numerical flux



.. _tut:poisson:2D:fig3:

.. figure:: eps/Poisson2D_Dvc_surf1.png
   :width: 400

   Examples on plots created by transforming the finite element field to a field on a uniform, structured 2D grid: a surface plot of the solution



.. figure:: eps/Poisson2D_Dvc_mesh1.png
   :width: 400

   Examples on plots created by transforming the finite element field to a field on a uniform, structured 2D grid: lifted mesh plot of the solution




It should be easy with the information above to transform a finite element
field over a uniform rectangular or box-shaped mesh to the corresponding
``BoxField`` object and perform Matlab-style
visualizations of the whole field or
the field over planes or along lines through the domain.
By the transformation to a regular grid we have some more flexibility
than what Viper offers. (It should be added that
comprehensive tools like
VisIt, MayaVi2, or ParaView also have the possibility for plotting fields
along lines and extracting planes in 3D geometries, though usually with
less degree of control compared to Gnuplot, Matlab, and Matplotlib.)

.. _tut:poisson:nD:

Parameterizing the Number of Space Dimensions
---------------------------------------------

.. index:: dimension-independent code


FEniCS makes it is easy to write a unified simulation code that can operate
in 1D, 2D, and 3D. We will conveniently make use of this feature in
forthcoming examples. The relevant technicalities are therefore explained
below.

Consider the simple problem

.. math::


        u''(x) = 2\hbox{ in }[0,1],\quad u(0)=0,\ u(1)=1,


with exact solution :math:`u(x)=x^2`. Our aim is to formulate and solve this
problem in a 2D and a 3D domain as well.
We may generalize the domain :math:`[0,1]` to a box of any size
in the :math:`y` and :math:`z` directions and pose homogeneous Neumann
conditions :math:`\partial u/\partial n = 0` at all additional boundaries
:math:`y=\mbox{const}` and :math:`z=\mbox{const}` to ensure that :math:`u` only varies with
:math:`x`. For example, let us choose
a unit hypercube as domain: :math:`\Omega = [0,1]^d`, where :math:`d` is the number
of space dimensions. The generalized $d$-dimensional Poisson problem
then reads

.. math::


          \begin{array}{rcll}
            \Delta u  &=  2 &\mbox{in } \Omega, \\
            u  &=  0 &\mbox{on } \Gamma_0,\\
            u  &=  1 &\mbox{on } \Gamma_1,\\
        {\partial u\over\partial n}  &=  0 &\mbox{on } \partial\Omega\backslash\left(
        \Gamma_0\cup\Gamma_1\right),
          \end{array}


where :math:`\Gamma_0` is the side of the hypercube where :math:`x=0`, and
where :math:`\Gamma_1` is the side where :math:`x=1`.

Implementing a PDE for any :math:`d` is no more
complicated than solving a problem with a specific number of dimensions.
The only non-trivial part of the code is actually to define the mesh.
We use the command line to provide user-input to the program. The first argument
can be the degree of the polynomial in the finite element basis functions.
Thereafter, we supply the
cell divisions in the various spatial directions. The number of
command-line arguments will then imply the number of space dimensions.
For example, writing ``3 10 3 4`` on the command line means that
we want to approximate :math:`u` by piecewise polynomials of degree 3,
and that the domain is a three-dimensional cube with :math:`10\times 3\times 4`
divisions in the :math:`x`, :math:`y`, and :math:`z` directions, respectively.
Each of the :math:`10\times 3\times 4 = 120` boxes will
be divided into six tetrahedra.
The Python code can be quite compact:

.. code-block:: python

        degree = int(sys.argv[1])
        divisions = [int(arg) for arg in sys.argv[2:]]
        d = len(divisions)
        domain_type = [UnitInterval, UnitSquare, UnitCube]
        mesh = domain_type[d-1](*divisions)
        V = FunctionSpace(mesh, 'CG', degree)

First note that although ``sys.argv[2:]`` holds the divisions of
the mesh, all elements of the list ``sys.argv[2:]`` are string objects,
so we need to explicitly convert each element to an integer.
The construction ``domain_type[d-1]`` will pick the right name of the
object used to define the domain and generate the mesh.
Moreover, the argument ``*divisions``
sends each component of the list ``divisions`` as a separate
argument. For example, in a 2D problem where ``divisions`` has
two elements, the statement

.. code-block:: python

        mesh = domain_type[d-1](*divisions)

is equivalent to

.. code-block:: python

        mesh = UnitSquare(divisions[0], divisions[1])


The next part of the program is to set up the boundary conditions.
Since the Neumann conditions have :math:`\partial u/\partial n=0` we can
omit the boundary integral from the weak form. We then only
need to take care of Dirichlet conditions at two sides:

.. code-block:: python

        tol = 1E-14   # tolerance for coordinate comparisons
        def Dirichlet_boundary0(x, on_boundary):
            return on_boundary and abs(x[0]) < tol

        def Dirichlet_boundary1(x, on_boundary):
            return on_boundary and abs(x[0] - 1) < tol

        bc0 = DirichletBC(V, Constant(0), Dirichlet_boundary0)
        bc1 = DirichletBC(V, Constant(1), Dirichlet_boundary1)
        bc = [bc0, bc1]

Note that this code is independent of the number of space dimensions.
So are the statements defining and solving
the variational problem:

.. code-block:: python

        v = TestFunction(V)
        u = TrialFunction(V)
        f = Constant(-2)
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx

        problem = VariationalProblem(a, L, bc)
        u = problem.solve()

The complete code is found in ``Poisson123D_DN1.py``.

Observe that if we actually want to test variations in one selected
space direction, parameterized by ``e``, we only need to
replace ``x[0]`` in the code by ``x[e]`` (!). The parameter
``e`` could be given as the second command-line argument.
This extension appears in the file ``Poisson123D_DN2.py``.
You can run a 3D problem with this code where :math:`u` varies in, e.g.,
:math:`z` direction and is approximated by, e.g., a 5-th degree polynomial.
For any legal input the numerical solution coincides with the
exact solution at the nodes (because the exact solution is a second-degree
polynomial).
