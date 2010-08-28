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

.. literalinclude:: demo.py
   :lines: 13

Some optimization options for the form compiler a prescribed:

.. literalinclude:: demo.py
   :lines: 15-20

The first line tells the form compiler to use C++ compiler optimizations when
compiling the generated code. The remainder is a dictionary of options which
will be passed to the form compiler. It lists the optimizations strategies
that we wish the form compiler to use when generating code.

A unit cube mesh, with 16 vertices in each direction, of tetrahedra is
created, and on this mesh a finite element space of continuous Lagrange
basis functions is defined:

.. literalinclude:: demo.py
   :lines: 23-24

.. index:: compiled subdomain

The portions of the boundary on which Dirichlet boundary conditions will be applied
are now defined:

.. literalinclude:: demo.py
   :lines: 26-28

The boundary subdomain ``left`` corresponds to the part of the boundary on
which :math:`x=0` and the boundary subdomain ``right`` corresponds to the
part of the boundary on which :math:`x=1`. Note that C++ syntax is used in
the ``compile_subdomains`` function since the function will be automatically
compiled into C++ code for efficiency. The variable ``on_boundary`` is true
for points on the boundary of a domain, and false otherwise.

.. index:: compiled expression

Functions for the Dirichlet boundary values are defined using compiled
expressions:

.. literalinclude:: demo.py
   :lines: 30-35

For the Expression ``r``, the Python dictionary named ``defaults`` is used
to automatically set values in the function string.

The boundary subdomains and the boundary condition expressions are collected
together in to ``DirichtletBC`` objects:

.. literalinclude:: demo.py
   :lines: 37-38

The function space ``V`` is included to indicate the functions to which
the boundary conditions should be applied.

Test and trial functions, and the most recent approximate solution :math:`u`
are define on the finite element space :math:`V`, and ``Constants`` are declared
for the body force (``B``) and traction (``T``) terms:

.. literalinclude:: demo.py
   :lines: 40-45

In place of ``Constant``, it is also possible to use ``as_vector``, e.g.
``B = as_vector( [0.0, -0.5, 0.0] )``. The advantage of ``Constant`` is that
it values can be changed without requiring re-generation and re-compilation
of C++ code. Using ``as_vector`` can eliminate some function calls during
assembly.

With the functions defined, the kinematic quantities involved in the model
are defined using UFL syntax:

.. literalinclude:: demo.py
   :lines: 47-54

The strain energy density and the total potential energy are now defined
using UFL syntax:

.. literalinclude:: demo.py
   :lines: 56-64

Like for the body force and traction vectors, ``Constant'' has been used
for the model parameters ``mu`` and ``lmbda`` to avoid re-generation of C++
code when changing model parameters. Note that ``lambda`` is a reserved
keyword in Python, hence the misspelling ``lmbda``.

Directional derivatives are now computed of :math:`\Pi` and :math:`L`
(see :eq:`first_variation` and :eq:`second_variation`):

.. literalinclude:: demo.py
   :lines: 66-70

The functional ``L`` and it Jacobian ``a``, together with the list of
Dirichet boundary condition ``[bcl, bcr]`` are then used to construct and
nonlinear variational problem which is then solved:

.. literalinclude:: demo.py
   :lines: 72-74

The argument ``nonlinear = True`` indicated that the problem is nonlinear and
that a Newton solver should be used. The earlier defined dictionary of from
compiler options is passed using ``form_compiler_parameters = ffc_options``

Finally, the solution is saved to a file named ``displacement.pvd`` in VTK
format, and the deformed mesh is plotted to the screen:

.. literalinclude:: demo.py
   :lines: 76-



Complete code
-------------

.. literalinclude:: demo.py
   :start-after: # Begin demo
