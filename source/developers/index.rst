.. Developer information.

.. _developers:

#####################
Developer information
#####################

This page contains information for FEniCS developers, including an
overview of the organization of the FEniCS project, how to use
Launchpad and Bazaar, and how to write code and documentation.

************
Organization
************

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. Each component is developed by one
or more authors. This means that each component can be developed at
its own pace. At the same time, we strive to make periodic and
coordinated releases of all components to ensure interoperability
between the components.

The FEniCS Project uses `Launchpad <http://www.launchpad.net>`_ as its
main development platform. An overview of all FEniCS projects on
Launchpad can be found `here
<https://launchpad.net/fenics-project>`__. :ref:`launchpad_pages` also
contains a collection of links to important Launchpad pages for the
various FEniCS components.

Each component of FEniCS is registered as a project on Launchpad and
each project is connected to a pair of teams. The first team consists
of developers and users interested in the development of the
project. This team is open to everyone and has a mailing list where
the development of the project is discussed. The second team (the core
team) manages the project and has write-access to the source code.

***************************
Instructions for developers
***************************

The following pages list important instructions for FEniCS developers.

.. toctree::
   :maxdepth: 1

   Taking part in the development <taking_part>
   Obtaining the source code <getting_code>
   Writing code <writing_code>
   contributing_code

.. _developers_checklist:

*******************************
Summary of developer guidelines
*******************************

Be kind to your fellow developers and follow these guidelines:

#. Subscribe to all relevant Launchpad mailinglists
   (:ref:`More details <developers_taking_part>`)

#. Write Code and documents according to the style guidelines
   (:ref:`More details <developers_writing_code>`)

#. Test before you push
   (:ref:`More details <before_committing>`)

#. Do not rewrite history
   (:ref:`More details <bzr_branch_workflow>`)

