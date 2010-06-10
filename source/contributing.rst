.. Notes on how to contribute to the FEniCS project.

.. _contributing:

############
Contributing
############

There are several ways to contribute to the FEniCS project.
The easiest way is to simply use the software and give feedback about the
experience on the appropriate mailing list.
However, since there are still plenty of things left to be done, we encourage
users to take more active part in developing software and writing
documentation.
Of course as an active developer it is also easier to influence the direction
an focus of the individual projects.
FEniCS uses `Launchpad <https://launchpad.net/>`_ for managing and hosting the
projects and `FEniCS development <http://www.fenics.org/wiki/Development>`_
contain more information about the organisation of projects on Launchpad.

*****************************
Helping out through Launchpad
*****************************

The core components of FEniCS all uses Launchpad for answers to questions from
users, specification tracking via blueprints, bug tracking and code hosting.

Answers
=======

New users of software of a given project often have questions regarding usage
and best practises.
They will (or at least they should) ask questions using the
`Launchpad Answers <https://help.launchpad.net/Answers>`_ functionality.
Of course the developers of the given project are in the best position to
answer these questions, but very often the number of developers is less than a
handful and since the developers are often active in more than one group or
project any help from you is greatly appreciated.

`DOLFIN questions <https://answers.launchpad.net/dolfin>`_ show the questions
and answers so far related to the DOLFIN project.

.. _contributing_blueprints:

Blueprints
==========

As a user you might sometimes be missing features in the projects which could
make your own code more efficient or simplify the way you do things in the
code.
The missing feature could even prevent you from solving the problem that you
want, or maybe you have found a better way of implementing a current feature.
FEniCS uses `Blueprints <https://help.launchpad.net/Blueprint>`_ to keep track
of the specifications for new features planned in the projects.
However, before adding a new blueprint it is a good idea to first discuss the
design on the relevant mailing list.

`DOLFIN blueprints <https://blueprints.launchpad.net/dolfin>`_ show the current
blueprints for the DOLFIN project.

Bugs
====

FEniCS software is under active development.
As a consequence bugs are likely to occur from time to time.
In the event that you encounter a bug, please file a bug report using the
Launchpad system for tracking `Bugs <https://help.launchpad.net/Bugs>`_.
Please include the most simple code example that exposes the bug.

`DOLFIN bugs <https://bugs.launchpad.net/dolfin>`_ show a list of currently
known bugs in the DOLFIN project.

Contributing code
=================

Finding bugs and drawing up blueprints is a good way to contribute to the
FEniCS project, however, at some point those bugs has to be fixed and the new
features should be implemented which requires contributing code.
The main repository, or branch, of for instance DOLFIN is owned by the
`FEniCS Core Team <https://launchpad.net/~dolfin-core>`_ which is a restricted
team.
Therefore, unless you are a member of the core team, you will not be able to
upload any code to the main repository directly.
Instead you will have to submit a :ref:`patch <contributing_patches>` or create
your own :ref:`branch <contributing_branches>` which will be merged into the
main branch by a member of the core team.

.. note::

    Membership of a core team is usually not granted until you have committed
    a series of high quality changesets over a longer period of time.

Coding style
------------

DOLFIN is written mainly in ``C++`` while projects like FFC and UFL are written
in ``Python``.
To streamline the source code and ease the job for maintainers that need to
read and edit large amounts of code, a style guide for developers is useful.
The FEniCS coding style of ``Python`` adheres to the
`PEP-8 <http://www.python.org/dev/peps/pep-0008/>`_ although it is not
strictly enforced.
:ref:`styleguides_cpp_coding_style` explains in more detail how ``C++`` code
for DOLFIN is supposed to look like.

Note that the above style guides are *guides* only, and they can be abandoned
in certain cases if readability demands it.

.. _contributing_license:

License
-------

FEniCS projects are released under the GNU GPL v2, GNU GPL v3 or
GNU LGPL v2.1 licenses.
Please see the relevant project page on Launchpad to find out which license
apply to a given project.
License and author information is put at the top of the files which you have
modified or added.
This information should be provided according to the following examples.

For ``C++`` (DOLFIN):

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

    __author__ = "Anders Logg (logg@simula.no)"
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

Usually, for simple bug fixes and minor changes, submitting a patch is simplest
method to get code uploaded to the main branch.
For instance, to create and submit a patch for DOLFIN the following procedure
should be applied:

#. Get the current development branch::

    $ bzr branch lp:dolfin

#. Modify the files.

#.  If your contribution consists of new files, add them to the correct
    location in the directory tree::

        $ bzr add <files>

    where ``<files>`` is the list of new files.
    You do not have to take any action for previously existing files 
    which have been modified. Do not add temporary or binary files.

#.  Update author and date information as described in the
    :ref:`contributing_license` section.

#.  Commit your changes::

        $ bzr ci -m "<description>"

    where ``<description>`` is a short description of what your patch
    accomplishes.

#.  Create a patch::

        $ bzr send -o dolfin-<identifier>-<date>.patch

    where ``<identifier>`` is a keyword that can be used to identify the patch
    as coming from you (your username, last name, first name, a nickname etc.)
    and ``<date>`` is today's date in the format ``yyyy-mm-dd``.
    This will create the patch ``dolfin-<identifier>-<date>.patch`` as a file
    in the current directory.

#.  Send the patch that you just created to the DOLFIN mailing list
    dolfin@lists.launchpad.net with a description of the patch.
    A developer with write access to the main repository will then review the
    code and merge it with the main branch (assuming that it passes the code
    review).

The procedure for creating a patch for other FEniCS projects is identical, only
the project name, ``dolfin``, is different.


.. _contributing_branches:

Branches
--------

When implementing new features, described in :ref:`contributing_blueprints`,
or fixing more involved bugs, development might take place over several days or
even months.
Instead of submitting a patch once coding is complete, it is a good idea to
create a new branch on Launchpad to let others follow the progress.
Using DOLFIN as an example, simply do the following:

#. Get the current development branch::

    $ bzr branch lp:dolfin

#.  Go to `Bazaar branches of DOLFIN <https://code.launchpad.net/dolfin>`_,
    click on ``Register a branch`` and follow the instructions.

#.  Start developing as usual and remember that regular commits makes it easier
    to follow the development.

#.  Push changesets to the new branch::

        $ bzr push lp:<path-to-branch-location>

    The first time you push to this location you should use the
    ``--use-existing-dir`` option.

#.  Once you have completed your work you should propose it for merging into
    the DOLFIN main branch.
    A developer with write access to the main repository will then review the
    code and merge it with the main branch (assuming that it passes the code
    review).

The procedure for using branches for other FEniCS projects is identical, only
the project name, ``dolfin``, is different.

*********************
Writing documentation
*********************

`FEniCS Documentation <https://launchpad.net/fenics-doc>`_ is a component of
the FEniCS project. It is therefore organised and maintained using the same
framework as other sub projects and you can :ref:`contribute <contributing>` in
the same way as to any other project.
FEniCS projects and particularly DOLFIN is under active development which means
that the documentation is constantly lagging behind. Any help to accommodate
this is greatly appreciated. In particular we need help to:

* keep existing documentation (:ref:`tutorial_index`, :ref:`demos_index`,
  :ref:`programmers_reference_index`) in sync with the development versions of
  FEniCS projects (syntax changes etc.);

* update the :ref:`programmers_reference_index` for new/deleted features in
  DOLFIN;

* add new demos showing new features or solving old problems using an
  alternative approach.

The documentation is generated by
`Sphinx <http://sphinx.pocoo.org/index.html>`_ and uses
`reStructuredText <http://docutils.sourceforge.net/rst.html>`_ (reST) as the
markup language.
`reStructuredText primer <http://sphinx.pocoo.org/rest.html>`_ and
`Sphinx Markup Constructs <http://sphinx.pocoo.org/markup/index.html>`_ are
good places to start.
The :ref:`styleguides_sphinx_coding_style` explains how the reST source files
should look like.

Documenting demos
=================

When adding a new demo to the documentation, or updating an existing one,
the below model should be followed:

* Summarise what features are demonstrated
* Problem and method description
* Explain how each step of the solution process is implemented (include code
  snippets if appropriate)
* Add complete source code files for download
* Link to relevant sections of the :ref:`programmers_reference_index` and to
  demos that shows alternative implementations (if any)
* Add keywords to the index
* All demos MUST be available in both ``C++`` and ``Python`` versions
* Have someone review the documentation

See the guide on how to
:ref:`document demos <styleguides_sphinx_documenting_demos>` for details on how
to implement each step and which files are needed.
Also refer to the :ref:`Poisson (C++ demo) <demos_cpp_pde_poisson>` and 
:ref:`Poisson (Python demo) <demos_python_pde_poisson>` as examples on how the
documentation should look like.

.. note::

    Currently, as we're migrating demos from DOLFIN, there are a lot of demos
    that needs documentation. Please see :ref:`demos_missing_demos` and
    consider lending a hand filling in the blanks.

Programmer's reference
======================

.. note::

    KBO: Figure out how this should be organised and what requirements we have.

    * Covers DOLFIN
    * Documenting classes/functions
    * Example code
    * Link to relevant demos for actual usage (use directive seealso?)
    * Integrate doc strings with PyDOLFIN
    * Code snippets (does not necessarily have to run)
    * Class index (*General Index* or *Global Module* index), use ``index``?
    * Put a link to a nicely documented class in the programmer's reference

Before committing your work
===========================

There are a few simple tests that should be run before committing your work on
the documentaion:

* Run the script ``build_docs`` in the top level directory to make sure that
  the documentation is successfully build without warnings
* Run the script ``verify_code_snippets.py`` in the top level directory to test
  that all code snippets in the demos are exact copies of the code available in
  the source code files.

.. note::

    KBO: There will probably be more things to put here later

Please fix any errors you might encounter running these scripts even if your
work did not introduce them or at least notify the mailing list
fenics@lists.launchpad.net in case you are unable to do so.

