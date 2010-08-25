.. Notes on how to contribute to the FEniCS Project.

.. _contributing:

############
Contributing
############

There are several ways to contribute to the FEniCS Project. The
easiest way is to simply use the software and give feedback about the
experience on the appropriate mailing list. However, since there are
still plenty of things left to be done, we encourage users to take
more active part in developing software and writing documentation.  As
an active developer, it is also easier to influence the direction and
focus of the project. The step from user to developer is just a patch
away.

*****************************
Helping out through Launchpad
*****************************

FEniCS uses `Launchpad <https://launchpad.net/>`_ for answering user
questions, mailing lists, specification tracking via blueprints, bug
tracking, and code hosting.

Answers
=======

New users often have questions regarding usage and best
practices. They will (or at least they should) ask questions using the
`Launchpad Answers <https://help.launchpad.net/Answers>`_
functionality. Naturally, the FEniCS developers are in a good position
to answer these questions, but as the developers often are busy
writing code (or answering other questions), any help from you is
greatly appreciated.

`FEniCS answers <https://answers.launchpad.net/fenics>`_ is the place
for general FEniCS questions and `DOLFIN answers
<https://answers.launchpad.net/dolfin>`_ is the place for questions
regarding the FEniCS user interface DOLFIN.

`Launchpad pages <launchpad_pages.html>`_ contains a collection of
links to the answer pages for all FEniCS components.

.. _contributing_blueprints:

Mailing lists
=============

Users are also welcome to take part in discussions on the FEniCS
mailing lists, but specific questions are better directed to the
relevant `Answer page <https://help.launchpad.net/Answers>`_ as
explained above.

We try to keep a friendly tone on our mailing lists, but sometimes
questions can only be settled by a heated debate. Don't be afraid to
step in and state your opinion.

The `FEniCS mailing list <mailto:fenics@lists.launchpad.net>`_ is the
place for general discussions regarding the FEniCS Project and the
`DOLFIN mailing list <mailto:dolfin@lists.launchpad.net>`_ is the
place for discussions regarding the FEniCS user interface DOLFIN.

`Launchpad pages <launchpad_pages.html>`_ contains a collection of
links to the mailing lists for all FEniCS components.

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

Bugs
====

FEniCS software is under active development. As a consequence, bugs
are likely to occur from time to time. In the event that you encounter
a bug, please file a bug report using the Launchpad system for
tracking `Bugs <https://help.launchpad.net/Bugs>`_.  Please include
the most simple code example that exposes the bug. Let's say that
again: *Please include the most simple code example that exposes the
bug*. No one will want to dig through your application code to find
the bug.

`FEniCS bugs <https://bugs.launchpad.net/fenics>`_ is the place for
general FEniCS bugs and `DOLFIN bugs
<https://bugs.launchpad.net/dolfin>`_ is the place for bugs in the
FEniCS user interface DOLFIN.

`Launchpad pages <launchpad_pages.html>`_ contains a collection of
links to the bug tracking pages for all FEniCS components.

Contributing code
=================

Finding bugs and drawing up blueprints is a good way to contribute to
the FEniCS Project. However, at some point those bugs have to be fixed
and new features need to be implemented which requires that someone
writes the code. The main repository, or branch, of for instance
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
