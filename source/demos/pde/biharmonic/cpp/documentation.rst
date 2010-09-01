.. Documentation for the biharmonic demo from DOLFIN.

.. _demos_cpp_pde_biharmonic:


Biharmonic equation
===================

.. include:: ../common.txt


Implementation
--------------

The implementation is split in two files, a form file containing the definition
of the variational forms expressed in UFL and the solver which is implemented
in a C++ file.

UFL form file
^^^^^^^^^^^^^

First we define the variational problem in UFL in the file called
:download:`Biharmonic.ufl`.

In the UFL file, the finite element space is defined:

.. code-block:: python

    # Elements
    element = FiniteElement("Lagrange", triangle, 2)

On the space ``element``, trial and test functions, and the source term
are defined:

.. code-block:: python

    # Trial and test functions
    u = TrialFunction(element)
    v = TestFunction(element)
    f = Coefficient(element)

Next, the outward unit normal to cell boundaries and a measure of the cell
size are defined. The average size of cells sharing a facet will be used
(``h_avg``).  The UFL syntax ``('+')`` and ``('-')`` restricts a function
to the ``('+')`` and ``('-')`` sides of a facet, respectively.
The penalty parameter ``alpha`` is made a ``Constant`` so that it can be changed
in the program without regenerating the code.

.. code-block:: python

    # Normal component, mesh size and right-hand side
    n  = element.cell().n
    h  = Constant(triangle)
    h_avg = (h('+') + h('-'))/2

    # Parameters
    alpha = Constant(triangle)

Finally the bilinear and linear forms are defined. Integrals over
internal facets are indicated by ``*dS``.

.. code-block:: python

    # Bilinear form
    a = inner(div(grad(u)), div(grad(v)))*dx \
      - inner(avg(div(grad(u))), jump(grad(v), n))*dS \
      - inner(jump(grad(u), n), avg(div(grad(v))))*dS \
      + alpha('+')/h_avg*inner(jump(grad(u), n), jump(grad(v),n))*dS

    # Linear form
    L = f*v*dx

Complete code
-------------

Complete UFL file
^^^^^^^^^^^^^^^^^

.. literalinclude:: Biharmonic.ufl
   :start-after: # Compile
   :language: python

Complete main file
^^^^^^^^^^^^^^^^^^

.. literalinclude:: main.cpp
   :start-after: // using
   :language: c++
