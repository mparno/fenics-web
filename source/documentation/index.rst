.. _documentation:

#############
Documentation
#############

*******************
The FEniCS Tutorial
*******************

A good starting point for new users is the :doc:`FEniCS tutorial
<tutorial/index>`. The tutorial will help you get quickly up and running
with solving differential equations in FEniCS. The tutorial focuses
exclusively on the FEniCS Python interface, since this is the simplest
approach to exploring FEniCS for beginners.

***************
The FEniCS Book
***************

:ref:`The FEniCS Book <book>`, *Automated Solution of Differential
Equations by the Finite Element Method*, is a comprehensive (700
pages) book documenting the mathematical methodology behind the FEniCS
Project and the software developed as part of the FEniCS Project. The
FEniCS Tutorial is included as the opening chapter of the FEniCS Book.

*****************
The FEniCS Manual
*****************

`The FEniCS Manual
<http://launchpad.net/fenics-book/trunk/final/+download/fenics-manual-2011-10-31.pdf>`__
is a 200-page excerpt from the FEniCS Book, including the FEniCS
Tutorial, an introduction to the finite element method and
documentation of DOLFIN and UFL.

*************
Release notes
*************

If you are updating your application code to a new FEniCS release,
make sure to check the :ref:`release notes <releases>` where you will
find detailed information about new features and interface changes.

.. _demos:

*****
Demos
*****

A simple way to build your first FEniCS application is to copy and
modify one of the existing demos from DOLFIN. The table below provides
links to demo documentation for different DOLFIN versions.

.. list-table::
    :widths: 100 50 50 50 175
    :header-rows: 0
    :class: center

    * - DOLFIN (Python demos)
      - 1.0.0
      - `1.0.beta2 <../documentation/dolfin/1.0.beta2/python/demo/index.html>`__
      - `1.0.beta <../documentation/dolfin/1.0.beta/python/demo/index.html>`__
      - `development version <../documentation/dolfin/dev/python/demo/index.html>`__
    * - DOLFIN (C++ demos)
      - 1.0.0
      - `1.0.beta2 <../documentation/dolfin/1.0.beta2/cpp/demo/index.html>`__
      - `1.0.beta <../documentation/dolfin/1.0.beta/cpp/demo/index.html>`__
      - `development version <../documentation/dolfin/dev/cpp/demo/index.html>`__

More demos can be found either in the ``demo`` directory of the DOLFIN
source tree or under one of the directories ``/usr/share/dolfin/demo``
(GNU/Linux),
``/Applications/FEniCS.app/Contents/Resources/share/dolfin/demo`` (Mac
OS X) or ``C:\FEniCS\share\dolfin\demo`` (Windows) if FEniCS was
installed from a binary package.

.. _programmers_references:

**********************
Programmer's reference
**********************

The FEniCS application programming interface (API) is documented
separately for each FEniCS component. The most important interfaces
are those of the C++/Python problem solving environment :ref:`DOLFIN
<about_components_dolfin>` and the form language :ref:`UFL
<about_components_ufl>`. The table below provides links to the API
documentation for different versions of DOLFIN, UFL and other FEniCS
components.

.. list-table::
    :widths: 100 50 50 50 175
    :header-rows: 0
    :class: center

    * - DOLFIN (Python API)
      - 1.0.0
      - `1.0.beta2 <../documentation/dolfin/1.0.beta2/python/genindex.html>`__
      - `1.0.beta <../documentation/dolfin/1.0.beta/python/genindex.html>`__
      - `development version <../documentation/dolfin/dev/python/genindex.html>`__
    * - DOLFIN (C++ API)
      - 1.0.0
      - `1.0.beta2 <../documentation/dolfin/1.0.beta2/cpp/genindex.html>`__
      - `1.0.beta <../documentation/dolfin/1.0.beta/cpp/genindex.html>`__
      - `development version <../documentation/dolfin/dev/cpp/genindex.html>`__
    * - UFL (API)
      - 1.0.0
      - `1.0-beta3 <../documentation/ufl/1.0-beta3/genindex.html>`__
      - `1.0-beta2 <../documentation/ufl/1.0-beta2/genindex.html>`__
      - `development version <../documentation/ufl/dev/genindex.html>`__

.. toctree::
   :hidden:

   tutorial/index
