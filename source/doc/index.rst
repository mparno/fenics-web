.. _documentation:

#############
Documentation
#############

********
Tutorial
********

A good starting point for new users is the :ref:`FEniCS tutorial
<doc_tutorial>`. The tutorial will help you get quickly up and running
with solving differential equations in FEniCS. The tutorial focuses
exclusively on the FEniCS Python interface, since this is the simplest
approach to exploring FEniCS for beginners.

**********************
Programmer's reference
**********************

The FEniCS application programming interface (API) is documented
separately for each FEniCS component. The most important interfaces
are those of the C++/Python problem solving environment :ref:`DOLFIN
<about_projects_dolfin>` and the form language :ref:`UFL
<about_projects_ufl>`. The table below provides links to the API
documentation for different versions of DOLFIN, UFL and other FEniCS
components.

.. list-table::
    :widths: 50 50 50 50 150
    :header-rows: 0
    :class: center

    * - DOLFIN
      - 1.0.0
      - 0.9.12
      - 0.9.11
      - development version
    * - UFL
      - 1.0.0
      - 0.9.2
      - 0.9.1
      - development version

.. toctree::
   :hidden:

   tutorial/index

*****
Demos
*****

A simple way to build your first FEniCS application is to copy and
modify one of the existing demos from DOLFIN. These demos can be found
either in the ``demo`` directory of the DOLFIN source tree or under
the directory ``/usr/share/doc/dolfin-doc/demo`` if you have installed
FEniCS from a binary package.

.. note::

   Where do the demos reside on Mac and Windows?
