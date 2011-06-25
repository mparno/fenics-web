.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)

.. _tut:poisson:nonlinear:

Nonlinear Problems
==================

Now we shall address how to solve nonlinear PDEs in FEniCS. Our
sample PDE for implementation is taken as a nonlinear Poisson equation:

.. math::


        -\nabla\cdot\left( q(u)\nabla u\right) = f\thinspace .


The coefficient :math:`q(u)` makes the equation nonlinear (unless :math:`q(u)`
is a constant).

To be able to easily verify our implementation,
we choose the domain, :math:`q(u)`, :math:`f`, and the boundary
conditions such that we have
a simple, exact solution :math:`u`. Let
:math:`\Omega` be the unit hypercube :math:`[0, 1]^d`
in :math:`d` dimensions, :math:`q(u)=(1+u)^m`, :math:`f=0`, :math:`u=0` for :math:`x_0=0`, :math:`u=1`
for :math:`x_0=1`, and :math:`\partial u/\partial n=0` at all other boundaries
:math:`x_i=0` and :math:`x_i=1`, :math:`i=1,\ldots,d-1`. The coordinates are now represented by
the symbols :math:`x_0,\ldots,x_{d-1}`. The exact solution is then

.. math::


        u(x_0,\ldots,x_d) = \left((2^{m+1}-1)x_0 + 1\right)^{1/(m+1)} - 1\thinspace .



The variational formulation of our model problem reads:
Find :math:`u \in V` such that

.. math::


          F(u; v) = 0 \quad \forall v \in \hat{V},


where

.. math::



        F(u; v) = \int_\Omega q(u)\nabla u\cdot \nabla v dx,


and

.. math::


            \hat{V} &= \{v \in H^1(\Omega) : v = 0 \mbox{ on } x_0=0\mbox{ and }x_0=1\}, \\
             V      &= \{v \in H^1(\Omega) : v = 0 \mbox{ on } x_0=0\mbox{ and } v = 1\mbox{ on }x_0=1\}\thinspace .


The discrete problem arises as usual by restricting :math:`V` and :math:`\hat V` to a
pair of discrete spaces. As usual, we omit any subscript on discrete
spaces and simply say :math:`V` and :math:`\hat V` are chosen finite dimensional
according to some mesh and element type.
The nonlinear problem then reads: Find :math:`u\in V` such that

.. math::


          F(u; v) = 0 \quad \forall v \in \hat{V},



with :math:`u = \sum_{j=1}^N U_j \phi_j`. Since :math:`F` is a nonlinear function
of :math:`u`, the variational statement gives rise to a system of
nonlinear algebraic equations.
From now on the interest is only in the discrete problem, and as mentioned
in the section :ref:`tut:poisson1:varform`,
we simply write :math:`u` instead of :math:`u_h` to get a closer notation between
the mathematics and the Python code. When the exact solution needs to
be distinguished, we denote it by :math:`u_{\mbox{e}}`.

FEniCS can be used in alternative ways for solving a nonlinear PDE
problem. We shall in the following subsections go through four
solution strategies:
1) a simple Picard-type iteration,
2) a Newton method at the algebraic level,
3) a Newton method at the PDE level, and
4) an automatic approach where FEniCS attacks the nonlinear variational
problem directly. The "black box" strategy 4) is definitely the
simplest one from a
programmer's point of view, but the others give more control of the
solution process for nonlinear equations (which also has some
pedagogical advantages).

.. _tut:nonlinear:Picard:

Picard Iteration
----------------

.. index:: Picard iteration


.. index:: successive substitutions


Picard iteration is an easy way of handling nonlinear PDEs: we simply
use a known, previous solution in the nonlinear terms so that these
terms become linear in the unknown :math:`u`. The strategy is also known as
the method of successive substitutions.
For our particular problem,
we use a known, previous solution in the coefficient :math:`q(u)`.
More precisely, given a solution :math:`u^k` from iteration :math:`k`, we seek a
new (hopefully improved) solution :math:`u^{k+1}` in iteration :math:`k+1` such
that :math:`u^{k+1}` solves the *linear problem*,

.. math::



        \nabla\cdot \left(q(u^k)\nabla u^{k+1}\right) = 0,\quad k=0,1,\ldots


The iterations require an initial guess :math:`u^0`.
The hope is that :math:`u^{k} \rightarrow u` as :math:`k\rightarrow\infty`, and that
:math:`u^{k+1}` is sufficiently close to the exact
solution :math:`u` of the discrete problem after just a few iterations.

We can easily formulate a variational problem for :math:`u^{k+1}` from
the last equation.
Equivalently, we can approximate :math:`q(u)` by :math:`q(u^k)` in
:math:`\int_\Omega q(u)\nabla u\cdot \nabla v dx`
to obtain the same linear variational problem.
In both cases, the problem consists of seeking
:math:`u^{k+1} \in V` such that

.. math::


          \tilde F(u^{k+1}; v) = 0 \quad \forall v \in \hat{V},\quad k=0,1,\ldots,


with

.. math::



        \tilde F(u^{k+1}; v) = \int_\Omega q(u^k)\nabla u^{k+1}\cdot \nabla v dx
        \thinspace .


Since this is a linear problem in the unknown :math:`u^{k+1}`, we can equivalently
use the formulation

.. math::


        a(u^{k+1},v) = L(v),


with

.. math::


        a(u,v) &= \int_\Omega q(u^k)\nabla u\cdot \nabla v dx
        \\
        L(v) &= 0\thinspace .



The iterations can be stopped when :math:`\epsilon\equiv ||u^{k+1}-u^k||
< \mbox{tol}`, where :math:`\mbox{tol}` is small, say :math:`10^{-5}`, or
when the number of iterations exceed some critical limit. The latter
case will pick up divergence of the method or unacceptable slow
convergence.

In the solution algorithm we only need to store :math:`u^k` and :math:`u^{k+1}`,
called ``uk`` and ``u`` in the code below.
The algorithm can then be expressed as follows:

.. code-block:: python

        def q(u):
            return (1+u)**m

        # Define variational problem
        v = TestFunction(V)
        u = TrialFunction(V)
        uk = interpolate(Expression('0.0'), V)  # previous (known) u
        a = inner(q(uk)*grad(u), grad(v))*dx
        f = Constant(0.0)
        L = f*v*dx

        # Picard iterations
        u = Function(V)     # new unknown function
        eps = 1.0           # error measure ||u-uk||
        tol = 1.0E-5        # tolerance
        iter = 0            # iteration counter
        maxiter = 25        # max no of iterations allowed
        while eps > tol and iter < maxiter:
            iter += 1
            problem = VariationalProblem(a, L, bc)
            u = problem.solve()
            diff = u.vector().array() - uk.vector().array()
            eps = numpy.linalg.norm(diff, ord=numpy.Inf)
            print 'Norm, iter=%d: %g' % (iter, eps)
            uk.assign(u)    # update for next iteration

We need to define the previous solution in the iterations, ``uk``,
as a finite element function so that ``uk`` can be updated with
``u`` at the end of the loop. We may create the initial
``Function`` `uk`
by interpolating
an ``Expression`` or a ``Constant``
to the same vector space as ``u`` lives in (``V``).

In the code above we demonstrate how to use
``numpy`` functionality to compute the norm of
the difference between the two most recent solutions. Here we apply
the maximum norm (:math:`\ell_\infty` norm) on the difference of the solution vectors
(``ord=1`` and ``ord=2`` give the :math:`\ell_1` and :math:`\ell_2` vector
norms -- other norms are possible for ``numpy`` arrays,
see ``pydoc numpy.linalg.norm``).

The file ``nlPoisson_Picard.py`` contains the complete code for
this problem. The implementation is :math:`d` dimensional, with mesh
construction and setting of Dirichlet conditions as explained in
the section :ref:`tut:poisson:nD`.
For a :math:`33\times 33` grid with :math:`m=2` we need 9 iterations for convergence
when the tolerance is :math:`10^{-5}`.

.. _tut:nonlinear:Newton:algebraic:

A Newton Method at the Algebraic Level
--------------------------------------

After having discretized our nonlinear PDE problem, we may
use Newton's method to solve the system of nonlinear algebraic equations.
From the continuous variational problem,
the discrete version results in a
system of equations for the unknown parameters :math:`U_1,\ldots, U_N`

.. math::



        F_i(U_1,\ldots,U_N) \equiv
        \sum_{j=1}^N
        \int_\Omega \left( q\left(\sum_{\ell=1}^NU_\ell\phi_\ell\right)
        \nabla \phi_j U_j\right)\cdot \nabla \hat\phi_i dx = 0,\quad i=1,\ldots,N\thinspace .


Newton's method for the system :math:`F_i(U_1,\ldots,U_j)=0`, :math:`i=1,\ldots,N`
can be formulated as

.. math::


        \sum_{j=1}^N
        {\partial \over\partial U_j} F_i(U_1^k,\ldots,U_N^k)\delta U_j
        &= -F_i(U_1^k,\ldots,U_N^k),\quad i=1,\ldots,N,\\
        U_j^{k+1} &= U_j^k + \omega\delta U_j,\quad j=1,\ldots,N,


where :math:`\omega\in [0,1]` is a relaxation parameter, and :math:`k` is
an iteration index. An initial guess :math:`u^0` must
be provided to start the algorithm.
The original Newton method has :math:`\omega=1`, but in problems where it is
difficult to obtain convergence,
so-called *under-relaxation* with :math:`\omega < 1` may help.

.. index:: under-relaxation


We need, in a program, to compute the Jacobian
matrix :math:`\partial F_i/\partial U_j`
and the right-hand side vector :math:`-F_i`.
Our present problem has :math:`F_i` given by above.
The derivative :math:`\partial F_i/\partial U_j` becomes

.. math::


        \int\limits_\Omega \left\lbrack
         q'(\sum_{\ell=1}^NU_\ell^k\phi_\ell)\phi_j
        \nabla (\sum_{j=1}^NU_j^k\phi_j)\cdot \nabla \hat\phi_i
        +
        q\left(\sum_{\ell=1}^NU_\ell^k\phi_\ell\right)
        \nabla \phi_j \cdot \nabla \hat\phi_i
        \right\rbrack
         dx\thinspace .



The following results were used to obtain the previous equation:

.. math::


        {\partial u\over\partial U_j} = {\partial\over\partial U_j}
        \sum_{j=1}^NU_j\phi_j = \phi_j,\quad {\partial\over\partial U_j}\nabla u = \nabla\phi_j,\quad {\partial\over\partial U_j}q(u) = q'(u)\phi_j\thinspace .


We can reformulate the Jacobian matrix
by introducing the short
notation :math:`u^k = \sum_{j=1}^NU_j^k\phi_j`:

.. math::


        {\partial F_i\over\partial U_j} =
        \int_\Omega \left\lbrack
        q'(u^k)\phi_j
        \nabla u^k \cdot \nabla \hat\phi_i
        +
        q(u^k)
        \nabla \phi_j \cdot \nabla \hat\phi_i
        \right\rbrack
         dx\thinspace .


In order to make FEniCS compute this matrix, we need to formulate a
corresponding variational problem. Looking at the
linear system of equations in Newton's method,

.. math::


        \sum_{j=1}^N {\partial F_i\over\partial U_j}\delta U_j = -F_i,\quad
        i=1,\ldots,N,


we can introduce :math:`v` as a general test function replacing :math:`\hat\phi_i`,
and we can identify the unknown
:math:`\delta u = \sum_{j=1}^N\delta U_j\phi_j`. From the linear system
we can now go "backwards" to construct the corresponding
discrete weak form

.. math::



        \int_\Omega \left\lbrack
        q'(u^k)\delta u
        \nabla u^k \cdot \nabla v
        +
        q(u^k)
        \nabla \delta u\cdot \nabla v
        \right\rbrack
         dx = - \int_\Omega q(u^k)
        \nabla u^k\cdot \nabla v dx\thinspace .


This equation fits the standard form
:math:`a(\delta u,v)=L(v)` with

.. math::


        a(\delta u,v) &=
        \int_\Omega \left\lbrack
        q'(u^k)\delta u
        \nabla u^k \cdot \nabla v
        +
        q(u^k)
        \nabla \delta u \cdot \nabla v
        \right\rbrack
         dx\\
        L(v) &= - \int_\Omega q(u^k)
        \nabla u^k\cdot \nabla v dx\thinspace .


Note the important feature in Newton's method
that the
previous solution :math:`u^k` replaces :math:`u`
in the formulas when computing the matrix
:math:`\partial F_i/\partial U_j` and vector :math:`F_i` for the linear system in
each Newton iteration.

We now turn to the implementation.
To obtain a good initial guess :math:`u^0`, we can solve a simplified, linear
problem, typically with :math:`q(u)=1`, which yields the standard Laplace
equation :math:`\Delta u^0 =0`. The recipe for solving this problem
appears in the sections :ref:`tut:poisson1:varform`,
:ref:`tut:poisson1:impl`, and :ref:`tut:poisson1:DN`.
The code for computing :math:`u^0` becomes as follows:

.. code-block:: python

        tol = 1E-14
        def left_boundary(x, on_boundary):
            return on_boundary and abs(x[0]) < tol

        def right_boundary(x, on_boundary):
            return on_boundary and abs(x[0]-1) < tol

        Gamma_0 = DirichletBC(V, Constant(0.0), left_boundary)
        Gamma_1 = DirichletBC(V, Constant(1.0), right_boundary)
        bc = [Gamma_0, Gamma_1]

        # Define variational problem for initial guess (q(u)=1, i.e., m=0)
        v = TestFunction(V)
        u = TrialFunction(V)
        a = inner(grad(u), grad(v))*dx
        f = Constant(0.0)
        L = f*v*dx
        A, b = assemble_system(a, L, bc_u)
        uk = Function(V)
        solve(A, uk.vector(), b)

Here, ``uk`` denotes the solution function for the previous
iteration, so that the solution
after each Newton iteration is ``u = uk + omega*du``.
Initially, ``uk`` is the initial guess we call :math:`u^0` in the mathematics.


The Dirichlet boundary conditions for the problem to be solved in each Newton
iteration are somewhat different than the conditions for :math:`u`.
Assuming that :math:`u^k` fulfills the
Dirichlet conditions for :math:`u`, :math:`\delta u` must be zero at the boundaries
where the Dirichlet conditions apply, in order for :math:`u^{k+1}=u^k + \omega\delta u` to fulfill
the right Dirichlet values. We therefore define an additional list of
Dirichlet boundary conditions objects for :math:`\delta u`:

.. code-block:: python

        Gamma_0_du = DirichletBC(V, Constant(0), LeftBoundary())
        Gamma_1_du = DirichletBC(V, Constant(0), RightBoundary())
        bc_du = [Gamma_0_du, Gamma_1_du]

The nonlinear coefficient and its derivative must be defined
before coding the weak form of the Newton system:

.. code-block:: python

        def q(u):
            return (1+u)**m

        def Dq(u):
            return m*(1+u)**(m-1)

        du = TrialFunction(V) # u = uk + omega*du
        a = inner(q(uk)*grad(du), grad(v))*dx + \
            inner(Dq(uk)*du*grad(uk), grad(v))*dx
        L = -inner(q(uk)*grad(uk), grad(v))*dx


The Newton iteration loop is very similar to the Picard iteration loop
in the section :ref:`tut:nonlinear:Picard`:

.. code-block:: python

        du = Function(V)
        u  = Function(V)  # u = uk + omega*du
        omega = 1.0       # relaxation parameter
        eps = 1.0
        tol = 1.0E-5
        iter = 0
        maxiter = 25
        while eps > tol and iter < maxiter:
            iter += 1
            A, b = assemble_system(a, L, bc_du)
            solve(A, du.vector(), b)
            eps = numpy.linalg.norm(du.vector().array(), ord=numpy.Inf)
            print 'Norm:', eps
            u.vector()[:] = uk.vector() + omega*du.vector()
            uk.assign(u)

There are other ways of implementing the
update of the solution as well:

.. code-block:: python

        u.assign(uk)  # u = uk
        u.vector().axpy(omega, du.vector())

        # or
        u.vector()[:] += omega*du.vector()

The ``axpy(a, y)`` operation adds a scalar ``a`` times a ``Vector``
``y`` to a ``Vector`` object.  It is usually a fast operation
calling up an optimized BLAS routine for the calculation.

Mesh construction for a $d$-dimensional problem with arbitrary degree of
the Lagrange elements can be done as
explained in the section :ref:`tut:poisson:nD`.
The complete program appears in the file ``nlPoisson_algNewton.py``.


.. _tut:nonlinear:Newton:pdelevel:

A Newton Method at the PDE Level
--------------------------------

Although Newton's method in PDE problems is normally formulated at the
linear algebra level, i.e., as a solution method for systems of nonlinear
algebraic equations, we can also formulate the method at the PDE level.
This approach yields a linearization of the PDEs before they are discretized.
FEniCS users will probably find this technique simpler to apply than
the more standard method of the section :ref:`tut:nonlinear:Newton:algebraic`.

Given an approximation to the solution field, :math:`u^k`, we seek a
perturbation :math:`\delta u` so that

.. math::


        u^{k+1} = u^k + \delta u


fulfills the nonlinear PDE.
However, the problem for :math:`\delta u` is still nonlinear and nothing is
gained. The idea is therefore to assume that :math:`\delta u` is sufficiently
small so that we can linearize the problem with respect to :math:`\delta u`.
Inserting :math:`u^{k+1}` in the PDE,
linearizing the :math:`q` term as

.. math::


        q(u^{k+1}) = q(u^k) + q'(u^k)\delta u + {\cal O}((\delta u)^2)
        \approx q(u^k) + q'(u^k)\delta u,


and dropping other nonlinear terms in :math:`\delta u`,
we get

.. math::


        \nabla\cdot\left( q(u^k)\nabla u^k\right) +
        \nabla\cdot\left( q(u^k)\nabla\delta u\right) +
        \nabla\cdot\left( q'(u^k)\delta u\nabla u^k\right) = 0\thinspace .


We may collect the terms with the unknown :math:`\delta u` on the left-hand side,

.. math::


        \nabla\cdot\left( q(u^k)\nabla\delta u\right) +
        \nabla\cdot\left( q'(u^k)\delta u\nabla u^k\right) =
        -\nabla\cdot\left( q(u^k)\nabla u^k\right),


The weak form of this PDE is derived by multiplying by a test function :math:`v`
and integrating over :math:`\Omega`, integrating the second-order derivatives
by parts:

.. math::


        \int_\Omega \left(
        q(u^k)\nabla\delta u\cdot \nabla v
        + q'(u^k)\delta u\nabla u^k\cdot \nabla v\right) dx
        = -\int_\Omega q(u^k)\nabla u^k\cdot \nabla v dx\thinspace .


The variational problem reads: Find :math:`\delta u\in V` such that
:math:`a(\delta u,v) = L(v)` for all :math:`v\in \hat V`, where

.. math::


        a(\delta u,v) &=
        \int_\Omega \left(
        q(u^k)\nabla\delta u\cdot \nabla v
        + q'(u^k)\delta u\nabla u^k\cdot \nabla v\right) dx,
        \\
        L(v) &= -
        \int_\Omega q(u^k)\nabla u^k\cdot \nabla v dx\thinspace .



The function spaces :math:`V` and :math:`\hat V`, being continuous or discrete,
are as in the
linear Poisson problem from the section :ref:`tut:poisson1:varform`.

We must provide some initial guess, e.g., the solution of the
PDE with :math:`q(u)=1`. The corresponding weak form :math:`a_0(u^0,v)=L_0(v)`
has

.. math::


        a_0(u,v)=\int_\Omega\nabla u\cdot \nabla v dx,\quad L(v)=0\thinspace .


Thereafter, we enter a loop and solve
:math:`a(\delta u,v)=L(v)` for :math:`\delta u` and compute a new approximation
:math:`u^{k+1} = u^k + \delta u`. Note that :math:`\delta u` is a correction, so if
:math:`u^0` satisfies the prescribed
Dirichlet conditions on some part :math:`\Gamma_D` of the boundary,
we must demand :math:`\delta u=0` on :math:`\Gamma_D`.

Looking at the equations just derived,
we see that the variational form is the same as for the Newton method
at the algebraic level in the section :ref:`tut:nonlinear:Newton:algebraic`. Since Newton's method at the
algebraic level required some "backward" construction of the
underlying weak forms, FEniCS users may prefer Newton's method at the
PDE level, which is more straightforward.  There is seemingly no need
for differentiations to derive a Jacobian matrix, but a mathematically
equivalent derivation is done when nonlinear terms are linearized
using the first two Taylor series terms and when products in the
perturbation :math:`\delta u` are neglected.

The implementation is identical to the one in
the section :ref:`tut:nonlinear:Newton:algebraic` and is found in
the file ``nlPoisson_pdeNewton.py`` (for the fun of it we use
a ``VariationalProblem`` object instead of assembling a matrix and
vector and calling ``solve``). The reader is encouraged to go
through this code to be convinced that the present method actually
ends up with the same program as needed for the Newton method at
the linear algebra level in the section :ref:`tut:nonlinear:Newton:algebraic`.


.. _tut:nonlinear:Newton:auto:

Solving the Nonlinear Variational Problem Directly
--------------------------------------------------

DOLFIN has a built-in Newton solver and is able to automate the
computation of nonlinear, stationary boundary-value problems.
The automation is demonstrated next. A nonlinear variational
problem
can be solved by

.. code-block:: python

        VariationalProblem(J, F, bc, nonlinear=True)

where ``F`` corresponds to the nonlinear form :math:`F(u;v)` and
``J`` is a form for the derivative of ``F``.

The appropriate ``F`` form
is straightforwardly defined (assuming ``q(u)`` is
coded as a Python function):

.. code-block:: python

        v = TestFunction(V)
        u = Function(V)  # the unknown
        F = inner(q(u)*grad(u), grad(v))*dx

Note here that ``u`` is a ``Function``, not a ``TrialFunction``.
We could, alternatively, define :math:`F(u;v)` directly in terms of
a trial function for :math:`u` and a test function for :math:`v`, and then
created the proper ``F`` by

.. code-block:: python

        v = TestFunction(V)
        u = TrialFunction(V)
        Fuv = inner(q(u)*grad(u), grad(v))*dx
        u = Function(V)  # previous guess
        F = action(Fuv, u)

The latter statement is equivalent to :math:`F(u=u_0; v)`, where :math:`u_0` is
an existing finite element function representing the most recently
computed approximation to the solution.


.. index:: Gateaux derivative


The derivative :math:`J` (``J``) of :math:`F` (``F``) is formally the
Gateaux derivative :math:`DF(u^k; \delta u, v)`
of :math:`F(u;v)` at :math:`u=u^k` in the direction of :math:`\delta u`.
Technically, this Gateaux derivative is derived by computing

.. math::


        \lim_{\epsilon\rightarrow 0}{d\over d\epsilon} F_i(u^k + \epsilon\delta u; v)
        \thinspace .


The :math:`\delta u` is now the trial function and :math:`u^k` is as usual the previous
approximation to the solution :math:`u`.
We start with

.. math::


        {d\over d\epsilon}\int_\Omega \nabla v\cdot\left( q(u^k + \epsilon\delta u)
        \nabla (u^k + \epsilon\delta u)\right) dx


and obtain

.. math::


        \int_\Omega \nabla v\cdot\left\lbrack
        q'(u^k + \epsilon\delta u)\delta u
        \nabla (u^k + \epsilon\delta u)
        +
        q(u^k + \epsilon\delta u)
        \nabla \delta u
        \right\rbrack dx,


which leads to

.. math::


        \int_\Omega \nabla v\cdot\left\lbrack
        q'(u^k)\delta u
        \nabla (u^k)
        +
        q(u^k)
        \nabla \delta u
        \right\rbrack dx,


as :math:`\epsilon\rightarrow 0`.
This last expression is the Gateaux derivative of :math:`F`. We may use :math:`J` or
:math:`a(\delta u, v)` for this derivative, the latter having the advantage
that we easily recognize the expression as a bilinear form. However, in
the forthcoming code examples ``J`` is used as variable name for
the Jacobian.
The specification of ``J`` goes as follows:

.. code-block:: python

        du = TrialFunction(V)
        J = inner(q(u)*grad(du), grad(v))*dx + \
            inner(Dq(u)*du*grad(u), grad(v))*dx

where ``u`` is a ``Function`` representing the most recent solution.


.. index:: derivative


The UFL language that we use to specify weak forms supports differentiation
of forms. This means that when ``F`` is given as above, we can simply
compute the Gateaux derivative by

.. code-block:: python

        J = derivative(F, u, du)

The differentiation is done symbolically so no numerical approximation
formulas are involved. The ``derivative`` function is obviously
very convenient in problems where differentiating ``F`` by hand
implies lengthy calculations.


.. index:: nonlinear variational problems


The solution of the nonlinear problem is now a question of two statements:

.. code-block:: python

        problem = VariationalProblem(J, F, bc, nonlinear=True)
        u = problem.solve(u)

The ``u`` we feed to ``problem.solve`` is filled with the solution and
returned, implying that the ``u`` on the left-hand side actually refers
to the same ``u`` as provided on the right-hand side.  Python has a
convention that all input data to a function or class method are
represented as arguments, while all output data are returned to the
calling code. Data used as both input and output, as in this case,
will then be arguments and returned. It is not necessary to have a
variable on the left-hand side, as the function object is modified
correctly anyway, but it is convention that we follow here.

The file ``nlPoisson_vp1.py`` contains the complete code where
``J`` is calculated manually, while ``nlPoisson_vp2.py`` is
a counterpart where ``J`` is computed by ``derivative(F, u, du)``.
The latter file represents clearly the most automated way of solving
the present nonlinear problem in FEniCS.
