.. Developer information.

.. _developer:

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
its own pace, but we strive to make periodic and coordinated releases
of all components to ensure interoperability between the components.

Initially, FEniCS consisted of just two components (DOLFIN and FIAT)
but over time, several new components have been added and FEniCS now
consists of more than 10 individual components. Some of these
components (such as FIAT and UFC) have matured and reached a more
stable state, while others are changing at a faster pace. Currently,
most development takes place in DOLFIN, the C++ and Python interface
of FEniCS.

***************
Using Launchpad
***************

The FEniCS Project uses `Launchpad <http://www.launchpad.net>`_ as the
main development platform.

Each component of FEniCS is registered as a project on Launchpad and
each project is connected to a pair of teams. The first team consists
of developers and users interested in the development of the
project. This team is open to everyone and has a mailing list where
the development of the project is discussed. The second team (the core
team) manages the project and has write-access to the source code.

The first thing a prospective developer needs to do is therefore to
register an account on Launchpad and add herself/himself to the
relevant project teams. Membership to core teams is only granted after
a developer has proven reliable by committing a significant number of
high quality contributions.

An overview of all FEniCS projects on Launchpad can be found `here
<https://launchpad.net/fenics-project>`_.  `Launchpad pages
<launchpad_pages.html>`_ also contains a collection of links to
important Launchpad pages for the various FEniCS components.

Below, we describe how Launchpad is used to handle user questions, bug
reports, blueprints, and code hosting.

User questions
==============

User questions are discussed on the :ref:`help_answers` pages.
Developers should make sure to join the relevant team for each component so
that they will be notified about new user questions.

Bug reports
===========

Bug reports are discussed on the :ref:`help_bugs` pages.
Developers should make sure to join the relevant team for each component so
that they will be notified about new bugs.

Blueprints
==========

FEniCS uses `blueprints <https://help.launchpad.net/Blueprint>`_ to
keep track of the specifications for new features. However, before
adding a new blueprint it is a good idea to first discuss the design
on the relevant mailing list.

`FEniCS blueprints <https://blueprints.launchpad.net/fenics>`_ is the
place for general FEniCS blueprints and `DOLFIN blueprints
<https://blueprints.launchpad.net/dolfin>`_ is the place for
blueprints regarding the FEniCS user interface DOLFIN.

`Launchpad pages <launchpad_pages.html>`_ contains a collection of
links to the blueprint pages for all FEniCS components.

Code hosting
============

FEniCS uses Launchpad for hosting code. Each FEniCS component has a
`Bazaar <http://bazaar.canonical.com/en/>`_ repository on Launchpad
that contains all source code (including the entire development
history) for the component. The repositories are readable for
everyone, but write access is only granted to the members of the core
teams.

Developers should subscribe to notifications of changes made to the
source code by visiting the repository on Launchpad and clicking the
subscribe button.

************
Using Bazaar
************

Here is a quick reference for `using Bazaar
<http://doc.bazaar-vcs.org/bzr.2.0/en/quick-reference/index.html>`_.
In addition, a few useful commands for Bazaar follow below.

To set your identity with Bazaar, type

.. code-block:: sh

    bzr whoami "My Name <myname@foo.com>"

To create a new branch:

.. code-block:: sh

    bzr branch <address-to-original-branch> [address-to-new-branch]

To commit a change:

.. code-block:: sh

    bzr commit

To push changes:

.. code-block:: sh

    bzr push <address-to-branch>

To pull changes:

.. code-block:: sh

    bzr pull <address-to-branch>

The current development version of each FEniCS component can be
obtained directly using a special shortcut for code hosted on
Launchpad:

.. code-block:: sh

    bzr branch lp:<project-name>

For instance, one may create a branch of the main DOLFIN repository by
typing

.. code-block:: sh

    bzr branch lp:dolfin

*****************
Contributing code
*****************

The main repository, or branch, for each FEniCS component is owned by
the core team of that component. Therefore, unless you are a member of
the core team, you will not be able to upload any code to the main
repository directly. Instead you will have to submit a :ref:`patch
<contributing_patches>`, or create your own :ref:`branch
<contributing_branches>`. If the code is accepted, the patch or branch
will be merged into the main branch by a member of the core team.


Style guides
============

To ease the job for maintainers that will need to read and understand
your code, read the :ref:`styleguides` that explain
how to format your code so that it matches the coding style used for
FEniCS.

Before committing your work
===========================

Before committing any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying.

.. _contributing_patches:

Creating a patch
================

For simple bug fixes and minor changes, submitting a patch is the
simplest method to get code uploaded to the main branch. For instance,
to create and submit a patch for DOLFIN, the following procedure
should be applied.

#. Get the current development branch::

    bzr branch lp:dolfin

#. Modify the files.

#. If your contribution consists of new files, add those to the
   repository::

    bzr add <files>

   where ``<files>`` is the list of new files. Do not add temporary or
   binary files. No action is necessary for previously existing files
   which have been modified.

#. Update the author and date information as described in the
   :ref:`license <contributing_license>` section.

#. Commit your changes::

    bzr commit -m "<description>"

   where ``<description>`` is a short description of what your patch
   accomplishes.

#. Create a patch::

    bzr send -o dolfin-<identifier>-<date>.patch

   where ``<identifier>`` is a keyword that can be used to identify
   the patch as coming from you (your username, last name, first name,
   a nickname etc.) and ``<date>`` is today's date in the format
   ``yyyy-mm-dd``. This will create the patch
   ``dolfin-<identifier>-<date>.patch`` as a file in the current
   directory.

#. Send the patch that you just created to the DOLFIN mailing list
   dolfin@lists.launchpad.net with a description of the patch. A
   developer with write access to the main repository will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for creating a patch for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

.. _contributing_branches:

Creating a branch
=================

When implementing new features, described in
:ref:`contributing_blueprints`, or fixing more involved bugs,
development might take place over several days or even months.
Instead of submitting a patch once coding is complete, it is a good
idea to create a new branch on Launchpad to let others follow the
progress.  Using DOLFIN as an example, simply do the following:

#. Get the current development branch::

    bzr branch lp:dolfin

#. Go to `DOLFIN code page on Launchpad
   <https://code.launchpad.net/dolfin>`_, click on **Register a
   branch**, and follow the instructions.

#. Start developing as usual and remember that regular commits make it
   easier to follow the development.

#. Push changesets to the new branch::

    bzr push lp:<path-to-branch-location>

   The first time you push to this location you should use the
   ``--use-existing-dir`` option.

#. Once you have completed your work, you should propose it for
   merging into the DOLFIN main branch (via the Launchpad system). A
   developer with write access to the main repository will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for using branches for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

*********************
Writing documentation
*********************

`FEniCS Documentation <https://launchpad.net/fenics-doc>`_ is a
component of the FEniCS Project. It is therefore organized and
maintained using the same framework as all other FEniCS components.
FEniCS and in particular DOLFIN are under active development, which
means that the documentation needs to be continuously updated. Any
help to accommodate this is greatly appreciated.

The documentation is generated by `Sphinx
<http://sphinx.pocoo.org/index.html>`_ and uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ (reST) as the markup
language.  Good starting points are `reStructuredText primer
<http://sphinx.pocoo.org/rest.html>`_ and `Sphinx Markup Constructs
<http://sphinx.pocoo.org/markup/index.html>`_.  The
:ref:`styleguides_sphinx_coding_style` explains what the reST source
files should look like.

.. _contributing_license:

***************
License
***************
