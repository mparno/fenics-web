.. Documentation for the mixed Poisson demo from DOLFIN.

.. _demos_cpp_pde_mixed-poisson:

Mixed formulation for Poisson's equation
========================================

.. include:: ../../../common/pde/mixed-poisson/mixed-poisson.txt


Implementation
--------------

The implementation is split in two files, a form file containing the definition
of the variational forms expressed in UFL and the solver which is implemented
in a C++ file.

Creating the form file
^^^^^^^^^^^^^^^^^^^^^^

First we define the variational problem in UFL which we save in the
file called :download:`MixedPoisson.ufl`.

We begin by defining the finite element spaces. We define two finite
element spaces :math:`\Sigma_h = BDM` and :math:`V_h = DG` separately,
before combining these into a mixed finite element space:

.. code-block:: python

    BDM = FiniteElement("BDM", "triangle", 1)
    DG  = FiniteElement("DG", "triangle", 0)
    W = BDM * DG

The first argument to ``FiniteElement`` specifies the type of finite
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

    (sigma, u) = TrialFunctions(W)
    (tau, v) = TestFunctions(W)

Further, we need to specify the source :math:`f` (a coefficient) that
will be used in the linear form of the variational problem. This
coefficient needs be defined on a finite element space, but none of
the above defined elements are quite appropriate. We therefore define
a separate finite element space for this coefficient.

.. code-block:: python

    CG = FiniteElement("CG", "triangle", 1)
    f = Coefficient(CG)

Finally, we define the bilinear and linear forms according to the equations:

.. code-block:: python

    a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
    L = - f*v*dx


Writing the solver
^^^^^^^^^^^^^^^^^^

The solver is implemented in the :download:`main.cpp` file.

At the top we include the DOLFIN header file and the header file containing the
variational forms for the Poisson equation.
For convenience we also include the DOLFIN namespace.

.. code-block:: c++

#include <dolfin.h>
#include "MixedPoisson.h"
