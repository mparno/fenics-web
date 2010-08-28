.. Documentation for the hyperelasticity demo from DOLFIN.

.. _demos_python_pde_hyperelasticity:

Hyperelasticity
===============

.. include:: ../../../common/pde/hyperelasticity/hyperelasticity.txt


Implementation
--------------

This demo is implemented in a single Python file, :download:`demo.py`, which
contains both the variational forms and the solver.

First, the ``dolfin`` module is imported:

.. code-block:: python

    from dolfin import *

Some optimization options for the form compiler a prescribed:

.. code-block:: python

    # Optimization options for the form
    parameters["form_compiler"]["cpp_optimize"] = True
    ffc_options = {"optimize": True, \
                   "eliminate_zeros": True, \
                   "precompute_basis_const": True, \
                   "precompute_ip_const": True}

The first line tells the form compiler to use C++ compiler optimizations when
compiling the generated code. The remainder is a dictionary of options which
will be passed to the form compiler. It lists the optimizations strategies
that we wish the form compiler to use when generating code.

A unit cube mesh, with 16 vertices in each direction, of tetrahedra is
created, and on this mesh a finite element space of continuous Lagrange
basis functions is defined:

.. code-block:: python

    # Create mesh and define function space
    mesh = UnitCube(16, 16, 16)
    V = VectorFunctionSpace(mesh, "Lagrange", 1)

.. index:: compiled subdomain

The portions of the boundary on which Dirichlet boundary conditions will be applied
are now defined:

.. code-block:: python

    # Mark boundary subdomians
    left, right = compile_subdomains(["(std::abs(x[0])       < DOLFIN_EPS) && on_boundary",
                                      "(std::abs(x[0] - 1.0) < DOLFIN_EPS) && on_boundary"])

The boundary subdomain ``left`` corresponds to the part of the boundary on
which :math:`x=0` and the boundary subdomain ``right`` corresponds to the
part of the boundary on which :math:`x=1`. Note that C++ syntax is used in
the ``compile_subdomains`` function since the function will be automatically
compiled into C++ code for efficiency. The variable ``on_boundary`` is true
for points on the boundary of a domain, and false otherwise.

.. index:: compiled expression

Functions for the Dirichlet boundary values are defined using compiled
expressions:

.. code-block:: python

    # Define Dirichlet boundary values (x = 0 or x = 1)
    c = Expression(("0.0", "0.0", "0.0"))
    r = Expression(("scale*0.0",
                    "scale*(y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta) - x[1])",
                    "scale*(z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta) - x[2])"),
                    defaults = dict(scale = 0.5, y0 = 0.5, z0 = 0.5, theta = pi/3))

For the Expression ``r``, the Python dictionary named ``defaults`` is used
to automatically set values in the function string.

The boundary subdomains and the boundary condition expressions are collected
together in to ``DirichtletBC`` objects:

.. code-block:: python

    bcl = DirichletBC(V, c, left)
    bcr = DirichletBC(V, r, right)

The function space ``V`` is included to indicate the functions to which
the boundary conditions should be applied.

Test and trial functions, and the most recent approximate solution :math:`u`
are define on the finite element space :math:`V`, and ``Constants`` are declared
for the body force (``B``) and traction (``T``) terms:

.. code-block:: python

    # Define functions
    v  = TestFunction(V)             # Test function
    du = TrialFunction(V)            # Incremental displacement
    u  = Function(V)                 # Displacement from previous iteration
    B  = Constant((0.0, -0.5, 0.0))  # Body force per unit mass
    T  = Constant((0.1,  0.0, 0.0))  # Traction force on the boundary

In place of ``Constant``, it is also possible to use ``as_vector``, e.g.
``B = as_vector( [0.0, -0.5, 0.0] )``. The advantage of ``Constant`` is that
it values can be changed without requiring re-generation and re-compilation
of C++ code. Using ``as_vector`` can eliminate some function calls during
assembly.

With the functions defined, the kinematic quantities involved in the model
are defined using UFL syntax:

.. code-block:: python

    # Kinematics
    I = Identity(V.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)

The strain energy density and the total potential energy are now defined
using UFL syntax:

.. code-block:: python

    # Elasticity parameters
    E, nu = 10.0, 0.3
    mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

    # Total potential energy
    Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

Like for the body force and traction vectors, ``Constant'' has been used
for the model parameters ``mu`` and ``lmbda`` to avoid re-generation of C++
code when changing model parameters. Note that ``lambda`` is a reserved
keyword in Python, hence the misspelling ``lmbda``.

Directional derivatives are now computed of :math:`\Pi` and :math:`L`
(see :eq:`first_variation` and :eq:`second_variation`):

.. code-block:: python

    # Compute first variation of Pi (directional derivative about u in the direction of v)
    L = derivative(Pi, u, v)

    # Compute Jacobian of L
    a = derivative(L, u, du)

The functional ``L`` and it Jacobian ``a``, together with the list of
Dirichet boundary condition ``[bcl, bcr]`` are then used to construct and
nonlinear variational problem which is then solved:

.. code-block:: python

    # Create nonlinear variational problem and solve
    problem = VariationalProblem(a, L, [bcl, bcr], nonlinear = True, form_compiler_parameters = ffc_options)
    problem.solve(u)

The argument ``nonlinear = True`` indicated that the problem is nonlinear and
that a Newton solver should be used. The earlier defined dictionary of from
compiler options is passed using ``form_compiler_parameters = ffc_options``

Finally, the solution is saved to a file named ``displacement.pvd`` in VTK
format, and the deformed mesh is plotted to the screen:

.. code-block:: python

    # Save solution in VTK format
    file = File("displacement.pvd");
    file << u;

    # Plot and hold solution
    plot(u, mode = "displacement", interactive = True)

Complete code
-------------

.. code-block:: python

    from dolfin import *

    # Optimization options for the form compiler
    parameters["form_compiler"]["cpp_optimize"] = True
    ffc_options = {"optimize": True, \
                   "eliminate_zeros": True, \
                   "precompute_basis_const": True, \
                   "precompute_ip_const": True}

    # Create mesh and define function space
    mesh = UnitCube(16, 16, 16)
    V = VectorFunctionSpace(mesh, "Lagrange", 1)

    # Mark boundary subdomians
    left, right = compile_subdomains(["(std::abs(x[0])       < DOLFIN_EPS) && on_boundary",
                                      "(std::abs(x[0] - 1.0) < DOLFIN_EPS) && on_boundary"])

    # Define Dirichlet boundary (x = 0 or x = 1)
    c = Expression(("0.0", "0.0", "0.0"))
    r = Expression(("scale*0.0",
                    "scale*(y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta) - x[1])",
                    "scale*(z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta) - x[2])"),
                    defaults = dict(scale = 0.5, y0 = 0.5, z0 = 0.5, theta = pi/3))

    bcl = DirichletBC(V, c, left)
    bcr = DirichletBC(V, r, right)

    # Define functions
    v  = TestFunction(V)             # Test function
    du = TrialFunction(V)            # Incremental displacement
    u  = Function(V)                 # Displacement from previous iteration
    B  = Constant((0.0, -0.5, 0.0))  # Body force per unit mass
    T  = Constant((0.1,  0.0, 0.0))  # Traction force on the boundary

    # Kinematics
    I = Identity(V.cell().d)    # Identity tensor
    F = I + grad(u)             # Deformation gradient
    C = F.T*F                   # Right Cauchy-Green tensor

    # Invariants of deformation tensors
    Ic = tr(C)
    J  = det(F)

    # Elasticity parameters
    E, nu = 10.0, 0.3
    mu, lmbda = Constant(E/(2*(1 + nu))), Constant(E*nu/((1 + nu)*(1 - 2*nu)))

    # Stored strain energy density (compressible neo-Hookean model)
    psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

    # Total potential energy
    Pi = psi*dx - dot(B, u)*dx - dot(T, u)*ds

    # Compute first variation of Pi (directional derivative about u in the direction of v)
    L = derivative(Pi, u, v)

    # Compute Jacobian of L
    a = derivative(L, u, du)

    # Create nonlinear variational problem and solve
    problem = VariationalProblem(a, L, [bcl, bcr], nonlinear = True, form_compiler_parameters = ffc_options)
    problem.solve(u)

    # Save solution in VTK format
    file = File("displacement.pvd");
    file << u;

    # Plot and hold solution
    plot(u, mode = "displacement", interactive = True)
