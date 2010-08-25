.. Developer information.

.. _developer:

#####################
Developer information
#####################


*****************************
Helping out through Launchpad
*****************************

FEniCS uses `Launchpad <https://launchpad.net/>`_ for answering user
questions, mailing lists, specification tracking via blueprints, bug
tracking, and code hosting.

There are several ways to contribute to the FEniCS Project. The
easiest way is to simply use the software and give feedback about the
experience on the appropriate mailing list. However, since there are
still plenty of things left to be done, we encourage users to take
more active part in developing software and writing documentation.  As
an active developer, it is also easier to influence the direction and
focus of the project. The step from user to developer is just a patch
away.

This page contains information for FEniCS developers.  FEniCS development
tools are hosted Launchpad. This page explains how to obtain source
code from `Launchpad <https://launchpad.net/>`_, subscribe to developer
mailing lists and take part in the development of FEniCS. An overview
of all FEniCS projects on Launchpad can be found
`here <https://launchpad.net/fenics-project>`_.

Blueprints
==========

As a user you might sometimes be missing features which could make
your own code more efficient or simplify the way you solve problems
with FEniCS. A missing feature could prevent you from solving the
problem that you want, or maybe you have found a better way of
implementing a current feature. FEniCS uses `Blueprints
<https://help.launchpad.net/Blueprint>`_ to keep track of the
specifications for new features. However, before adding a new
blueprint it is a good idea to first discuss the design on the
relevant mailing list.

`FEniCS blueprints <https://blueprints.launchpad.net/fenics>`_ is the
place for general FEniCS blueprints and `DOLFIN blueprints
<https://blueprints.launchpad.net/dolfin>`_ is the place for
blueprints regarding the FEniCS user interface DOLFIN.

`Launchpad pages <launchpad_pages.html>`_ contains a collection of
links to the blueprint pages for all FEniCS components.

**********************************
Organization of FEniCS development
**********************************

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. Each project is registered as a project
on Launchpad and each project is connected to a pair of teams. The first
project team consists of developers and users interested in the development
of the project. This team is open to everyone and has a mailing list where
the development of the project is discussed. The second team (the core team)
manages the project and has write-access to the source code.

********************************
Getting the development versions
********************************

The source code of all FEniCS projects can be obtained directly from Launchpad
using Bazaar, the version control system used by FEniCS. To get the source
code for a FEniCS project, simply do

.. code-block:: sh

    bzr branch lp:project-name

For example, to obtain the source code for DOLFIN, run the command

.. code-block:: sh

    bzr branch lp:dolfin

*************
Mailing lists
*************

To subscribe to the mailing list for a project, you must first join the project
team. The project teams are open to everyone. To join a project, visit the
project team page and click "Join the team" in the top right corner. You may
then click the link "Subscribe to mailing list" in the bottom left corner
of the project team page.

A link to the project team page can be found in the table below.

Note that mailing lists do not automatically receive bug reports and
notifications of branch changes (commit logs). See below for instructions
on how to subscribe to bug reports and branch changes.

Note also that core team members are automatically subscribed to bug mail and
branch notifications, but not to mailing lists so developers must manually
subscribe to the project mailing list.

**************************
Subscribing to bug reports
**************************

To subscribe to bug reports for a project, visit the main Launchpad page of
the project and click the link "Subscribe to bug mail" in the top right corner.

*****************************
Subscribing to branch changes
*****************************

To subscribe to branch changes for a project, visit the main Launchpad
page of the project and click the link "Branches" in the top menu. Then
click on the branch (lp::foo) in the list and finally "Subscribe yourself"
in the top right corner.

************
Using Bazaar
************

A quick reference for using Bazaar can be found
`here <http://doc.bazaar-vcs.org/bzr.2.0/en/quick-reference/index.html>`_.
To set your identity with Bazaar, type

.. code-block:: sh

    bzr whoami "My Name <myname@foo.com>"

To create a new branch (similar to hg clone):

.. code-block:: sh

    bzr branch address-to-original-branch [address-to-new-branch]

To branch a project hosted by Launchpad:

.. code-block:: sh

    bzr branch lp:project-name

To commit a change:

.. code-block:: sh

    bzr commit

To push changes:

.. code-block:: sh

    bzr push [address-to-branch]

To pull changes:

.. code-block:: sh

    bzr pull [address-to-branch]



 The main repository, or branch, of for instance
DOLFIN is owned by the `DOLFIN Core Team
<https://launchpad.net/~dolfin-core>`_ which is a restricted team.
Therefore, unless you are a member of the core team, you will not be
able to upload any code to the main repository directly. Instead you
will have to submit a :ref:`patch <contributing_patches>` or create
your own :ref:`branch <contributing_branches>` which will be merged
into the main branch by a member of the core team (if accepted).

.. note::

    Membership of a core team is usually not granted until you have
    committed a series of high quality patches/changesets over a
    longer period of time.

Coding style
------------

The FEniCS user interface DOLFIN is written mainly in C++ while other
FEniCS components like FFC, FIAT, and UFL are written in Python. To
streamline the source code and ease the job for maintainers that need
to read and edit large amounts of code, a style guide for developers
is useful. The FEniCS coding style for Python code adheres to the
`PEP-8 style guide <http://www.python.org/dev/peps/pep-0008/>`_
although it is not strictly enforced.
:ref:`styleguides_cpp_coding_style` explains in more detail the
preferred coding style for C++ code.

Note that the above style guides are *guides* only, and they can be
abandoned in certain cases if readability demands it.

.. _contributing_license:

License
-------

FEniCS components are released under the GNU GPL v2, GNU GPL v3, or
GNU LGPL v2.1 licenses. Please see the relevant component page on
Launchpad to find out which license applies to a given
component. License and author information is put at the top of the
files which you have modified or added. This information should be
provided according to the following examples.

For C++ (DOLFIN):

.. code-block:: c++

    // Copyright (C) 2007-2009 Anders Logg.
    // Licensed under the GNU LGPL Version 2.1.
    //
    // Modified by Garth N. Wells, 2007-2008.
    // Modified by Ola Skavhaug, 2008.
    //
    // First added:  2007-01-17
    // Last changed: 2009-06-22

For ``Python``:

.. code-block:: python

    __author__ = "Anders Logg <logg@simula.no>"
    __date__ = "2007-02-05"
    __copyright__ = "Copyright (C) 2007-2010 " + __author__
    __license__  = "GNU GPL version 3 or any later version"

    # Modified by Kristian B. Oelgaard, 2010.
    # Modified by Dag Lindbo, 2008.
    # Modified by Garth N. Wells, 2009.
    # Last changed: 2010-01-24

.. _contributing_patches:

Patches
-------

Usually, for simple bug fixes and minor changes, submitting a patch is
the simplest method to get code uploaded to the main branch. For
instance, to create and submit a patch for DOLFIN, the following
procedure should be applied:

#. Get the current development branch::

    bzr branch lp:dolfin

#. Modify the files.

#. If your contribution consists of new files, add them to the correct
   location in the directory tree::

    bzr add <files>

   where ``<files>`` is the list of new files. You do not have to take
   any action for previously existing files which have been
   modified. Do not add temporary or binary files.

#. Update the author and date information as described in the
   :ref:`contributing_license` section.

#. Commit your changes::

    bzr ci -m "<description>"

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

Branches
--------

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




Contributing code
=================

Finding bugs and drawing up blueprints is a good way to contribute to
the FEniCS Project. However, at some point those bugs have to be fixed
and new features need to be implemented which requires that someone
writes the code. For information about how to contribute code to the
FEniCS project, see the `<Developer page <developer>`_.

*********************
Writing documentation
*********************

`FEniCS Documentation <https://launchpad.net/fenics-doc>`_ is a
component of the FEniCS Project. It is therefore organized and
maintained using the same framework as all other FEniCS components and
you can :ref:`contribute <contributing>` in the same way as to any
other component. FEniCS components and in particular DOLFIN is under
active development, which means that the documentation needs to be
continuously updated. Any help to accommodate this is greatly
appreciated. In particular, we need help to:

The documentation is generated by `Sphinx
<http://sphinx.pocoo.org/index.html>`_ and uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ (reST) as the markup
language.  Good starting points are `reStructuredText primer
<http://sphinx.pocoo.org/rest.html>`_ and `Sphinx Markup Constructs
<http://sphinx.pocoo.org/markup/index.html>`_.  The
:ref:`styleguides_sphinx_coding_style` explains what the reST source
files should look like.

Programmer's reference
======================

The programmer's reference covers the C++ and Python interfaces to
FEniCS with emphasis on the DOLFIN library. The documentation of a
class/function should in general follow the below structure:

* One line which summarizes the funtionality of the class/function
* *Arguments*, a description of arguments
* *Returns*, a description of return values
* *Example*, a short code snippets that illustrate basic usage. The code does not
  have to be a stand-alone program.
* *See also*, links to demos which use the particular feature which is
  documented.

See the guide on how to :ref:`document the FEniCS interface
<styleguides_sphinx_documenting_interface>` for details on how to
implement each step and which files are needed.  Also refer to the
:ref:`Mesh class C++ documentation
<programmers_reference_cpp_mesh_Mesh>` and :ref:`Mesh class Python
documentation <programmers_reference_python_mesh_Mesh>` for good
examples of what the documentation should look like.

Documenting demos
=================

When adding a new demo to the documentation, or updating an existing one,
the below model should be followed:

* Summarize what features are demonstrated
* Problem and method description
* Explain how each step of the solution process is implemented (include code
  snippets if appropriate)
* Add complete source code files for download
* Link to relevant sections of the :ref:`programmers_reference_index` and to
  demos that show alternative implementations (if any)
* Add keywords to the index
* Make the demo available in both C++ and Python versions (this is important!)
* Have someone review the documentation

See the guide on how to :ref:`document demos
<styleguides_sphinx_documenting_demos>` for details on how to
implement each step and which files are needed.  Also refer to the
:ref:`Poisson C++ demo <demos_cpp_pde_poisson>` and :ref:`Poisson
Python demo <demos_python_pde_poisson>` for good examples.

.. note::

    Currently, as we're migrating demos from the DOLFIN source tree
    into this documentation, there are many demos that need
    documentation. Please see :ref:`demos_missing_demos` and consider
    lending a hand to in the blanks.

Before committing your work
===========================

There are a few simple tests that should be run before committing your work on
the documentaion:

* Run the script ``test/verify_demo_code_snippets.py`` to test that all code
  snippets in the demos are exact copies of the code available in the source
  code files.
* Run ``make all`` in the top level directory to make sure that
  the documentation is successfully build without warnings

Please fix any errors you might encounter running these scripts even
if your work did not introduce them or at least notify the mailing
list fenics@lists.launchpad.net in case you are unable to do so.
