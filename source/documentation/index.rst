.. _documentation:

#############
Documentation
#############

********
Tutorial
********

A good starting point for new users is the :doc:`FEniCS tutorial
<tutorial/index>`. The tutorial will help you get quickly up and running
with solving differential equations in FEniCS. The tutorial focuses
exclusively on the FEniCS Python interface, since this is the simplest
approach to exploring FEniCS for beginners.

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
    :widths: 100 50 50 50 50 175
    :header-rows: 0
    :class: center

    * - DOLFIN (Python API)
      - 1.0.0
      - `1.0.beta <http://www.fenicsproject.org/new/doc/dolfin/1.0.beta/python/genindex.html>`__
      - 0.9.12
      - 0.9.11
      - `development version <http://www.fenicsproject.org/new/doc/dolfin/dev/python/genindex.html>`__
    * - DOLFIN (C++ API)
      - 1.0.0
      - `1.0.beta <http://www.fenicsproject.org/new/doc/dolfin/1.0.beta/cpp/genindex.html>`__
      - 0.9.12
      - 0.9.11
      - `development version <http://www.fenicsproject.org/new/doc/dolfin/dev/cpp/genindex.html>`__
    * - UFL (API)
      - 1.0.0
      - `1.0-beta2 <http://www.fenicsproject.org/new/doc/ufl/1.0-beta2/genindex.html>`__
      - 0.9.2
      - 0.9.1
      - `development version <http://www.fenicsproject.org/new/doc/ufl/dev/genindex.html>`__

.. toctree::
   :hidden:

   tutorial/index

.. note::
    Figure out how to set relative paths rather than absolute paths as above.

.. note::
    Replace index.html by genindex.html in above links

.. note::
    Make page title more verbose for each link, like "DOLFIN C++ API version x.y.z"

*****
Demos
*****

A simple way to build your first FEniCS application is to copy and
modify one of the existing demos from DOLFIN. The table below provides
links to demo documentation for different DOLFIN versions.

.. list-table::
    :widths: 100 50 50 50 50 175
    :header-rows: 0
    :class: center

    * - DOLFIN (Python demos)
      - 1.0.0
      - `1.0.beta <http://www.fenicsproject.org/new/doc/dolfin/1.0.beta/python/demo>`__
      - 0.9.12
      - 0.9.11
      - `development version <http://www.fenicsproject.org/new/doc/dolfin/dev/python/demo>`__
    * - DOLFIN (C++ demos)
      - 1.0.0
      - `1.0.beta <http://www.fenicsproject.org/new/doc/dolfin/1.0.beta/cpp/demo>`__
      - 0.9.12
      - 0.9.11
      - `development version <http://www.fenicsproject.org/new/doc/dolfin/dev/cpp/demo>`__

More demos can be found either in the ``demo`` directory of the DOLFIN
source tree or under one of the directories ``/usr/share/dolfin/demo``
(GNU/Linux),
``/Applications/FEniCS.app/Contents/Resources/share/dolfin/demo`` (Mac
OS X) or ``C:\FEniCS\share\dolfin\demo`` (Windows) if FEniCS was
installed from a binary package.
