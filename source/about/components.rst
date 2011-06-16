.. _about_components:

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. A short overview of these components
is given below. Further information can be found in the :ref:`FEniCS
Book <book>` and on the :ref:`Launchpad page <launchpad_pages>` for
each individual component.

.. _about_components_core:

###############
Core components
###############

.. _about_components_dolfin:

******
DOLFIN
******

DOLFIN is a C++/Python library that functions as the main user
interface of FEniCS. A large part of the functionality of FEniCS is
implemented as part of DOLFIN. It provides a problem solving
environment for models based on partial differential equations and
implements core parts of the functionality of FEniCS, including data
structures and algorithms for computational meshes and finite element
assembly. To provide a simple and consistent user interface, DOLFIN
wraps the functionality of other FEniCS components and external
software, and handles the communication between these components.

Maintainers (in alpabetical order)
==================================

Johan Hake, Anders Logg, Garth N. Wells

Authors (past and present in chronological order)
=================================================

Anders Logg, Johan Hoffman, Garth N. Wells, Johan Jansson, Ola
Skavhaug, Kent-Andre Mardal, Martin Sandve Alnes, Johan Hake, Niclas
Jansson, Johannes Ring, Kristian Oelgaard, Marie Rognes

For a full list of contributors, refer to the ``AUTHORS`` file found
in the DOLFIN source tree.

.. note::
    Is this a good way to give credit? Who is missing? Is the order correct?

.. _about_components_ffc:

***
FFC
***

One of the key features of FEniCS is automated code generation for the
general and efficient solution of finite element variational
problems. FFC (FEniCS Form Compiler) is a compiler for variational
forms. It generates efficient low-level C++ code (UFC) from a
high-level mathematical description (UFL) of a finite element
variational problem.

.. image:: images/ufl-ffc-ufc.png
    :align: center

Maintainers (in alpabetical order)
==================================

Anders Logg, Kristian Oelgaard, Marie Rognes, Garth N. Wells

Authors (past and present in chronological order)
=================================================

Anders Logg, Marie Rognes, Kristian Oelgaard, Garth N. Wells

.. _about_components_fiat:

****
FIAT
****

FIAT (FInite element Automatic Tabulator) supports generation of
arbitrary order instances of the Lagrange elements on lines,
triangles, and tetrahedra. It is also capable of generating arbitrary
order instances of Jacobi-type quadrature rules on the same element
shapes. Further, H(div) and H(curl) conforming finite element spaces
such as the families of Raviart-Thomas, Brezzi-Douglas-Marini and
Nedelec are supported on triangles and tetrahedra. Upcoming versions
will also support Hermite and nonconforming elements.

Maintainers (in alpabetical order)
==================================

Robert C. Kirby

Authors (past and present in chronological order)
=================================================

Robert C. Kirby, Anders Logg, Marie Rognes

.. _about_components_instant:

*******
Instant
*******

Instant is a Python module that allows for instant inlining of C and
C++ code in Python. It is a small Python module built on top of SWIG
and Distutils. Instant is used by FFC and DOLFIN for just-in-time
(JIT) compilation of variational forms and expressions.

Maintainers (in alpabetical order)
==================================

???

Authors (past and present in chronological order)
=================================================

???

.. _about_components_ufc:

***
UFC
***

UFC (Unified Form-assembly Code) is a unified framework for finite
element assembly. More precisely, it defines a fixed interface for
communicating low level routines (functions) for evaluating and
assembling finite element variational forms. The UFC interface
consists of a single header file ufc.h that specifies a C++ interface
that must be implemented by code that complies with the UFC
specification.

Maintainers (in alpabetical order)
==================================

Anders Logg, Martin Sandve Alnes, Garth N. Wells

Authors (past and present in chronological order)
=================================================

Anders Logg, Martin Sandve Alnes, Kent-Andre Mardal, Ola Skavhaug,
Hans Petter Langtangen, Garth N. Wells

.. _about_components_ufl:

***
UFL
***

UFL (Unified Form Language) is a domain specific language for
declaration of finite element discretizations of variational
forms. More precisely, it defines a flexible interface for choosing
finite element spaces and defining expressions for weak forms in a
notation close to mathematical notation.

Maintainers (in alpabetical order)
==================================

Martin Sandve Alnes

Authors (past and present in chronological order)
=================================================

Martin Sandve Alnes, Anders Logg, Garth N. Wells

.. _about_components_additional:

#####################
Additional components
#####################

.. _about_components_ascot:

*****
ASCoT
*****

.. _about_components_dorsal:

Maintainers (in alpabetical order)
==================================

Marie Rognes

Authors (past and present in chronological order)
=================================================

Marie Rognes

******
Dorsal
******

Dorsal is a set of simple scripts to build components of the FEniCS
Project (as well as their dependencies) for various platforms.

.. _about_components_syfi:

Maintainers (in alpabetical order)
==================================

Harish Narayanan

Authors (past and present in chronological order)
=================================================

Harish Narayanan

********
SyFi/SFC
********

Maintainers (in alpabetical order)
==================================

Kent-Andre Mardal, Martin Sandve Alnes

Authors (past and present in chronological order)
=================================================

Kent-Andre Mardal, Martin Sandve Alnes

.. _about_components_viper:

*****
Viper
*****

Viper is a minimalistic scientific plotter and run-time visualization
module based on VTK. If installed, Viper provides built-in plotting
for DOLFIN. [`read more <https://launchpad.net/fenics-viper>`__]

Maintainers (in alpabetical order)
==================================

Ola Skavhaug

Authors (past and present in chronological order)
=================================================

Ola Skavhaug

.. _about_components_applications:

############
Applications
############

`FEniCS Apps <https://answers.launchpad.net/fenics-group>`__ is a
collection of applications/solvers based on FEniCS. If you have
developed a FEniCS-based application that you think qualifies to be on
this list, contact the `FEniCS Apps maintainers
<https://launchpad.net/~fenics-apps-core>`__.

* CBC.RANS
* CBC.Block
* `CBC.Solve <https://launchpad.net/cbc.solve>`__,
  a collection of Python-based PDE solvers
* DOLFWAVE
* DiffSim
* FEniCS Plasticity
* `GenFoo <https://launchpad.net/genfoo>`__,
  a generalized Fokker-Planck solver
* PUM Compiler
* PUM solver
* rheagen
* TriTetMesh
* `Unicorn <https://launchpad.net/unicorn>`__,
  a unified continuum mechanics solver

.. note::
    Which should be included among core components? Should we include SyFi?
    Viper? FErari? Those are strictly speaking not necessary to run an application.

.. note::
    Which applications are missing?

.. note::
    Should ASCoT be a component or an application?

.. note::
    Everyone should review the text presented for each component.

.. note::
    When we're happy with the information listed here, we should
    update the corresponding text on Launchpad.

.. note::
    Add some pretty pictures.
