.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)

.. _tut:linsys:

Controlling the Solution of Linear Systems
==========================================



Several linear algebra packages, referred to as
linear algebra *backends*, can be used in FEniCS to solve
linear systems:
PETSc, uBLAS, Epetra (Trilinos), or MTL4.
Which backend to apply can be controlled by setting

.. code-block:: python

        parameters['linear algebra backend'] = backendname

where ``backendname`` is a string, either ``'PETSc'``,
``'uBLAS'``, ``'Epetra'``, or ``'MTL4'``.
These backends offer high-quality implementations of both iterative
and direct solvers for linear systems of equations.


.. index:: down-casting matrices and vectors


The backend determines the specific data structures that are
used in the ``Matrix`` and ``Vector`` classes. For example,
with the PETSc backend, ``Matrix`` encapsulates a
PETSc matrix storage structure, and
``Vector`` encapsulates a PETSc vector storage structure.
Sometimes one wants to perform operations directly on (say) the
underlying PETSc objects. These can be fetched by

.. code-block:: python

        A_PETSc = down_cast(A).mat()
        b_PETSc = down_cast(b).vec()
        U_PETSc = down_cast(u.vector()).vec()

Here, ``u`` is a ``Function``, ``A`` is a ``Matrix``,
and ``b`` is a ``Vector``.
The same syntax applies if we want to fetch
the underlying Epetra, uBLAS, or MTL4 matrices and vectors. The section :ref:`tut:Epetra` provides an example on working directly with
Epetra objects.


Let us explain how one can choose between direct and iterative solvers.
We have seen that there
are two ways of solving linear systems, either we call the ``solve()``
method in a ``VariationalProblem`` object or we call the ``solve(A, U, b)``
function with the assembled coefficient matrix ``A``,
right-hand side vector ``b``, and
solution vector ``U``.

Variational Problem Objects
---------------------------

In case we use a ``VariationalProblem`` object, named ``problem``,
it has a ``parameters`` object that behaves like a Python dictionary,
and we can use this object to choose between a direct or iterative
solver:

.. code-block:: python

        problem.parameters['linear_solver'] = 'direct'
        # or
        problem.parameters['linear_solver'] = 'iterative'

Another parameter ``'symmetric'`` can be set to ``True``
if the coefficient matrix is symmetric so that
a method exploiting symmetry can be utilized.
For example, the default iterative solver is GMRES, but when solving
a Poisson equation, the iterative solution process will be more
efficient by setting the ``'symmetry'`` parameter
so that a Conjugate Gradient method is applied.

Having chosen an iterative solver, we can invoke
a submenu ``'krylov_solver'``
in the
``parameters`` object for setting various parameters for
the iterative solver (GMRES or Conjugate Gradients, depending on
whether the matrix is symmetric or not):

.. code-block:: python

        itsolver = problem.parameters['krylov_solver'] # short form
        itsolver['absolute_tolerance'] = 1E-10
        itsolver['relative_tolerance'] = 1E-6
        itsolver['divergence_limit'] = 1000.0
        itsolver['gmres_restart'] = 50
        itsolver['monitor_convergence'] = True
        itsolver['report'] = True

Here, ``'divergence_limit'``
governs the maximum allowable number of iterations,

the ``'gmres_restart'`` parameter tells how many iterations GMRES performs before
it restarts,
``'monitor_convergence'`` prints detailed information about the
development of the residual of a solver,
``'report'`` governs whether a one-line report about the solution
method and the number of iterations
is written on the screen or not. The absolute and relative tolerances
enter (usually residual-based) stopping criteria, which are dependent on
the implementation of the underlying iterative solver in the actual backend.



When direct solver is chosen, there is similarly a submenu
``'lu_solver'`` to set parameters, but here only the ``'report'``
parameter is available (since direct solvers very soldom have any
adjustable parameters). For nonlinear problems there is also
submenu ``'newton_solver'`` where tolerances, maximum iterations, and
so on, for a the Newton solver in ``VariationalProblem`` can be set.


.. index:: info


A complete list of all parameters and their default values
is printed to the screen by

.. code-block:: python

        info(problem.parameters, True)


Solve Function
--------------


.. index:: solve(A, x, b)


For the ``solve(A, U, b)`` approach, a 4th argument to ``solve``
determines the type of method:

  * ``'lu'`` for a sparse direct (LU decomposition) method,

  * ``'cg'`` for the Conjugate Gradient (CG) method, which is

applicable if ``A`` is symmetric and positive definite,
  * ``'gmres'`` for

the GMRES iterative method, which is applicable when ``A`` is nonsymmetric,
  * ``'bicgstab'`` for

the BiCGStab iterative method, which is applicatble when ``A`` is
nonsymmetric.

The default solver is ``'lu'``.

Good performance of an iterative method requires preconditioning of
the linear system. The 5th argument to ``solve`` determines the
preconditioner:

  * ``'none'`` for no preconditioning.

  * ``'jacobi'`` for the simple Jacobi (diagonal) preconditioner,

  * ``'sor'`` for SOR preconditioning,

  * ``'ilu'`` for incomplete LU factorization (ILU) preconditioning,

  * ``'icc'`` for incomplete Cholesky factorization preconditioning

(requires ``A`` to be symmetric and positive definite),
  * ``'amg_hypre'`` for algebraic multigrid (AMG) preconditioning

with the Hypre package (if available),
  * ``'mag_ml'`` for algebraic multigrid (AMG) preconditioning

with the ML package from Trilinos (if available),
  * ``'default_pc'`` for a default preconditioner, which depends

on the linear algebra backend (``'ilu'`` for PETSc).
.. , a domain decomposition with ILU subdomain solver for Epetra).

If the 5th argument is not provided, ``'ilu'`` is taken as the default
preconditioner.

Here are some sample calls to ``solve`` demonstrating the choice
of solvers and preconditioners:

.. code-block:: python

        solve(A, u.vector(), b)         # 'lu' is default solver
        solve(A, u.vector(), b, 'cg')   # CG with ILU prec.
        solve(A, u.vector(), b, 'gmres', 'amg_ml')  # GMRES with ML prec.



Setting the Start Vector
------------------------

.. index:: start vector for linear solvers


The choice of start vector for the iterations in a linear solver is often
important. With the ``solve(A, U, b)`` function the start vector
is the vector we feed in for the solution. A start vector
with random numbers in the interval :math:`[-1,1]` can be computed as

.. code-block:: python

        n = u.vector().array().size
        u.vector()[:] = numpy.random.uniform(-1, 1, n)
        solve(A, u.vector(), b, 'cg', 'ilu')

Or if a ``VariationalProblem`` object is used, its ``solve``
method may take an optional ``u`` function as argument (which we
can fill with the right values):

.. code-block:: python

        problem = VariationalProblem(a, L, bc)
        n = u.vector().array().size
        u.vector()[:] = numpy.random.uniform(-1, 1, n)
        u = problem.solve(u)


The program ``Poisson2D_DN_laprm.py`` demonstrates the various control mechanisms for
steering linear solvers as described above.

.. _tut:Epetra:

Using a Backend-Specific Solver
-------------------------------

.. ../../../la/trilinos/python/demo.py

Sometimes one wants to implement tailored solution algorithms, using
special features of the underlying numerical packages.
Here is an example where we create an ML preconditioned Conjugate
Gradient solver by programming with Trilinos-specific objects directly.
Given a linear system
:math:`AU=b`, represented by a ``Matrix`` object ``A``,
and two ``Vector`` objects ``U`` and ``b`` in a
Python program, the purpose is to
set up a solver using the Aztec Conjugate Gradient method from
Trilinos' Aztec library and combine that solver with the
algebraic multigrid preconditioner ML
from the ML library in Trilinos. Since the various parts of
Trilinos are mirrored in Python through the PyTrilinos package,
we can operate directly
on Trilinos-specific objects.

.. code-block:: python

        try:
            from PyTrilinos import Epetra, AztecOO, TriUtils, ML
        except:
            print '''You Need to have PyTrilinos with'
        Epetra, AztecOO, TriUtils and ML installed
        for this demo to run'''
            exit()

        from dolfin import *

        if not has_la_backend('Epetra'):
            print 'Warning: Dolfin is not compiled with Trilinos'
            exit()

        parameters['linear_algebra_backend'] = 'Epetra'

        # create matrix A and vector b in the usual way
        # u is a Function

        # Fetch underlying Epetra objects
        A_epetra = down_cast(A).mat()
        b_epetra = down_cast(b).vec()
        U_epetra = down_cast(u.vector()).vec()

        # Sets up the parameters for ML using a python dictionary
        ML_param = {"max levels"        : 3,
                    "output"            : 10,
                    "smoother: type"    : "ML symmetric Gauss-Seidel",
                    "aggregation: type" : "Uncoupled",
                    "ML validate parameter list" : False
        }

        # Create the preconditioner
        prec = ML.MultiLevelPreconditioner(A_epetra, False)
        prec.SetParameterList(ML_param)
        prec.ComputePreconditioner()

        # Create solver and solve system
        solver = AztecOO.AztecOO(A_epetra, U_epetra, b_epetra)
        solver.SetPrecOperator(prec)
        solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
        solver.SetAztecOption(AztecOO.AZ_output, 16)
        solver.Iterate(MaxIters=1550, Tolerance=1e-5)

        plot(u)
