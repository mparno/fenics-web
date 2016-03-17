.. _about_components:

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. A short overview of these components
is given below. Further information can be found in the :ref:`FEniCS
Book <book>` and on the `Bitbucket page
<https://bitbucket.org/fenics-project>`__.

.. _about_components_core:

###############
Core components
###############

.. _about_components_dolfin:

********************************************************
`DOLFIN <https://bitbucket.org/fenics-project/dolfin>`__
********************************************************

DOLFIN is a C++/Python library that functions as the main user
interface of FEniCS. A large part of the functionality of FEniCS is
implemented as part of DOLFIN. It provides a problem solving
environment for models based on partial differential equations and
implements core parts of the functionality of FEniCS, including data
structures and algorithms for computational meshes and finite element
assembly. To provide a simple and consistent user interface, DOLFIN
wraps the functionality of other FEniCS components and external
software, and handles the communication between these components.


Maintainers (in alphabetical order)
===================================

Johan Hake, Anders Logg, Chris Richardson, Garth N. Wells


Authors (past and present in chronological order)
=================================================

Anders Logg, Johan Hoffman, Garth N. Wells, Johan Jansson, Ola
Skavhaug, Kent-Andre Mardal, Martin Sandve Alnes, Johan Hake, Niclas
Jansson, Johannes Ring, Kristian B. Ølgaard, Marie Rognes

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/dolfin/raw/master/AUTHORS>`__ in
the DOLFIN source tree.

.. _about_components_ffc:


**************************************************
`FFC <https://bitbucket.org/fenics-project/ffc>`__
**************************************************

One of the key features of FEniCS is automated code generation for the
general and efficient solution of finite element variational
problems. FFC (FEniCS Form Compiler) is a compiler for variational
forms. It generates efficient low-level C++ code (UFC) from a
high-level mathematical description (UFL) of a finite element
variational problem.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <img src="../_images/ufl-ffc-ufc.png" class="img-responsive" alt="UFL-FFC-UFC">
    </div>
  </div>

Maintainers (in alphabetical order)
===================================

Martin Sandve Alnes, Anders Logg, Kristian B. Ølgaard, Marie Rognes,
Garth N. Wells


Authors (past and present in chronological order)
=================================================

Martin Sandve Alnes, Anders Logg, Marie Rognes, Kristian Oelgaard,
Garth N. Wells

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/ffc/raw/master/AUTHORS>`__ in
the FFC source tree.

.. _about_components_fiat:


****************************************************
`FIAT <https://bitbucket.org/fenics-project/fiat>`__
****************************************************

FIAT (FInite element Automatic Tabulator) supports generation of
arbitrary order instances of the Lagrange elements on lines,
triangles, and tetrahedra. It is also capable of generating arbitrary
order instances of Jacobi-type quadrature rules on the same element
shapes. Further, H(div) and H(curl) conforming finite element spaces
such as the families of Raviart-Thomas, Brezzi-Douglas-Marini and
Nedelec are supported on triangles and tetrahedra. Upcoming versions
will also support Hermite and nonconforming elements.

Maintainers (in alphabetical order)
===================================

Robert C. Kirby


Authors (past and present in chronological order)
=================================================

Robert C. Kirby, Anders Logg, Marie Rognes

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/fiat/raw/master/AUTHORS>`__ in
the FIAT source tree.


.. _about_components_instant:

**********************************************************
`Instant <https://bitbucket.org/fenics-project/instant>`__
**********************************************************

Instant is a Python module that allows for instant inlining of C and
C++ code in Python. It is a small Python module built on top of SWIG
and Distutils. Instant is used by FFC and DOLFIN for just-in-time
(JIT) compilation of variational forms and expressions.

Maintainers (in alphabetical order)
===================================

Martin Sandve Alnes


Authors (past and present in chronological order)
=================================================

.. note::
   Add list of Instant authors.

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/instant/raw/master/AUTHORS>`__
in the Instant source tree.

.. _about_components_ufc:

*************************************************************
`UFC <https://bitbucket.org/fenics-project/ufc-deprecated>`__
*************************************************************

UFC (Unified Form-assembly Code) is a unified framework for finite
element assembly. More precisely, it defines a fixed interface for
communicating low level routines (functions) for evaluating and
assembling finite element variational forms. The UFC interface
consists of a single header file ufc.h that specifies a C++ interface
that must be implemented by code that complies with the UFC
specification.

Maintainers (in alphabetical order)
===================================

Anders Logg, Martin Sandve Alnes, Garth N. Wells

Authors (past and present in chronological order)
=================================================

Anders Logg, Martin Sandve Alnes, Kent-Andre Mardal, Ola Skavhaug,
Hans Petter Langtangen, Garth N. Wells

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/ufc-deprecated/raw/master/AUTHORS>`__ in
the UFC source tree.

.. _about_components_ufl:

**************************************************
`UFL <https://bitbucket.org/fenics-project/ufl>`__
**************************************************

UFL (Unified Form Language) is a domain specific language for
declaration of finite element discretizations of variational
forms. More precisely, it defines a flexible interface for choosing
finite element spaces and defining expressions for weak forms in a
notation close to mathematical notation.

Maintainers (in alphabetical order)
===================================

Martin Sandve Alnes

Authors (past and present in chronological order)
=================================================

Martin Sandve Alnes, Anders Logg, Garth N. Wells, Kristian B. Ølgaard,
Marie E. Rognes

For a full list of contributors, refer to the file `AUTHORS
<https://bitbucket.org/fenics-project/ufl/raw/master/AUTHORS>`__ in
the UFL source tree.

.. _about_components_additional:
