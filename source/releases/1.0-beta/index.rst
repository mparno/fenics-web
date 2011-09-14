.. _release_1_0_beta:

#################################
Release notes for FEniCS 1.0-beta
#################################

The release of :ref:`FEniCS 1.0-beta <featured_1_0_beta>` is an
important step towards the release of FEniCS 1.0. A number of
important improvements were made with the release of 1.0-beta compared
to the previous release (0.9.11). These changes mostly involve `bug
fixes and code cleanups
<http://bazaar.launchpad.net/~dolfin-core/dolfin/main/view/head:/ChangeLog>`__.

Users will find that the interface for solving variational problems
has changed in 1.0-beta. We appreciate that this is an inconvenience
to users but believe that this is a necessary change. The interface
changes are detailed below.

Removal of the VariationalProblem class
=======================================

The change involves removal of the class ``VariationalProblem`` which
has been replaced by two new classes ``LinearVariationalProblem`` and
``NonlinearVariationalProblem``. In addition, two new classes
``LinearVariationalSolver`` and ``NonlinearVariationalSolver`` have
been added.

The following example illustrates the changes made for a linear
variational problem which, using the old interface, was implemented as
follows:

.. code-block:: python

    problem = VariationalProblem(a, L, bcs)
    solver.parameters["linear_solver"] = "gmres"
    solver.parameters["preconditioner"] = "ilu"
    u = problem.solve()

Using FEniCS 1.0-beta, the above syntax has changed to:

.. code-block:: python

    u = Function(V)
    problem = LinearVariationalProblem(a, L, u, bcs=bcs)
    solver = LinearVariationalSolver(problem)
    solver.parameters["linear_solver"] = "gmres"
    solver.parameters["preconditioner"] = "ilu"
    solver.solve()

Similarly, the syntax for a nonlinear variational problem now reads:

.. code-block:: python

    u = Function(V)
    problem = NonlinearVariationalProblem(F, u, bcs=bcs, J=J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters["linear_solver"] = "gmres"
    solver.parameters["preconditioner"] = "ilu"
    solver.solve()

Here, ``J`` is an optional argument that specifies the Jacobian of the
nonlinear form ``F``. If the Jacobian is not specified, it will be
automatically computed.

The new solve( ) function
=========================

As a shortcut for the above examples using
``LinearVariationalProblem`` and ``NonlinearVariationalProblem``, one
may instead use the new ``solve()`` function. The new solve function can
be used to solve linear systems as before:

.. code-block:: python

    solve(A, x, b)

However, it can also be used to solve linear and nonlinear variational
problems. A linear variational problem can be solved as follows:

.. code-block:: python

    u = Function(V)
    solve(a == L, u, bcs=bcs)

Similarly, a nonlinear variational problem can be solved as follows:

.. code-block:: python

    u = Function(V)
    solve(F == 0, u, bcs=bcs)

Solver and preconditioner arguments
===================================

The solver and preconditioner parameters to the ``solve()`` function
for linear systems have undergone a cleanup. For example, to solve a
linear system using GMRES with an algebraic multigrid preconditioner,
use

.. code-block:: python

    solve(A, x, b, "gmres", "amg")

The following tables list possible values for the solver and
preconditioner parameters:

Solver parameters
-----------------

+------------------+----------------------------------------------+
| Usage            | Method                                       |
+==================+==============================================+
| ``"lu"``         | LU factorization                             |
+------------------+----------------------------------------------+
| ``"cholesky"``   | Cholesky factorization                       |
+------------------+----------------------------------------------+
| ``"cg"``         | Conjugate gradient method                    |
+------------------+----------------------------------------------+
| ``"gmres"``      | Generalized minimal residual method          |
+------------------+----------------------------------------------+
| ``"bicgstab"``   | Biconjugate gradient stabilized method       |
+------------------+----------------------------------------------+
| ``"minres"``     | Minimal residual method                      |
+------------------+----------------------------------------------+
| ``"tfqmr"``      | Transpose-free quasi-minimal residual method |
+------------------+----------------------------------------------+
| ``"richardson"`` | Richardson method                            |
+------------------+----------------------------------------------+

Preconditioner parameters
-------------------------

+------------------------+----------------------------------------------+
| Usage                  | Preconditioner                               |
+========================+==============================================+
| ``"none"``             | No preconditioner                            |
+------------------------+----------------------------------------------+
| ``"ilu"``              | Incomplete LU factorization                  |
+------------------------+----------------------------------------------+
| ``"icc"``              | Incomplete Cholesky factorization            |
+------------------------+----------------------------------------------+
| ``"jacobi"``           | Jacobi iteration                             |
+------------------------+----------------------------------------------+
| ``"bjacobi"``          | Block Jacobi iteration                       |
+------------------------+----------------------------------------------+
| ``"sor"``              | Successive over-relaxation                   |
+------------------------+----------------------------------------------+
| ``"amg"``              | Algebraic multigrid (BoomerAMG or ML)        |
+------------------------+----------------------------------------------+
| ``"additive_schwarz"`` | Additive Schwarz                             |
+------------------------+----------------------------------------------+
| ``"hypre_amg"``        | Hypre algebraic multigrid (BoomerAMG)        |
+------------------------+----------------------------------------------+
| ``"hypre_euclid"``     | Hypre parallel incomplete LU factorization   |
+------------------------+----------------------------------------------+
| ``"hypre_parasails"``  | Hypre parallel sparse approximate inverse    |
+------------------------+----------------------------------------------+
| ``"ml_amg"``           | ML algebraic multigrid                       |
+------------------------+----------------------------------------------+

Default arguments required for Expressions
==========================================

The class ``Expression`` now requires default values for variables
used to define the expression. Thus, the following example:

.. code-block:: python

    f = Expression("sin(c*t)")
    ...
    f.c = 1.0
    f.t = 0.0

must be replaced by

.. code-block:: python

    f = Expression("sin(c*t)", c=1.0, t=1.0)
    ...
    f.c = 1.0
    f.t = 0.0
