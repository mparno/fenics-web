.. _developers_contributing_code:

*****************
Contributing code
*****************

The main branch for each FEniCS component is owned by the core team of
that component. Therefore, unless you are a member of the core team,
you will not be able to upload any code to the main branch
directly. Instead you will have to create your own :ref:`branch
<contributing_branches>` or submit a :ref:`patch
<contributing_patches>`. If the code is accepted, the patch or branch
will be merged into the main branch by a member of the core team.

If you are not familiar with `Bazaar <http://bazaar.canonical.com>`__,
the distributed revision control system used for all FEniCS
components, a good starting point is `Bazaar in five minutes
<http://doc.bazaar.canonical.com/latest/en/mini-tutorial>`__ or the
`Bazaar user guide
<http://doc.bazaar.canonical.com/latest/en/user-guide/index.html>`__.

.. _before_committing:

Before committing your work
===========================

Before committing any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying. If you are
submitting code for a new feature, it is expected that the new feature
is accompanied by a suitable set of unit tests. You should also review
the :ref:`recommended coding style <developers_writing_code>`.

.. _copyright_consent:

Copyright and licensing consent
===============================

Before your contribution can be accepted into FEniCS, you must sign a
copyright consent form. Ideally, both you and your employer should
sign a form. The forms can be found `here
<http://fenicsproject.org/pub/copyright/forms/>`__. After you have
signed the form, send it by regular mail to

  | Johannes Ring
  | Simula Research Laboratory
  | PO Box 134
  | 1325 Lysaker
  | Norway

Copies of signed consent forms are archived
`here <http://fenicsproject.org/pub/copyright/authors>`__
and `here <http://fenicsproject.org/pub/copyright/institutions>`__.

.. _contributing_branches:

Creating a branch
=================

A simple way to submit your changes is to create a branch on Launchpad
and submit a merge request. This allows the maintainers to pull the
source code from your branch, review it and then push to the main
development branch. Using DOLFIN as an example, simply do the
following:

#. Get the current development branch::

    bzr branch lp:dolfin

#. Go to `DOLFIN code page on Launchpad
   <https://code.launchpad.net/dolfin>`_, click on **Register a
   branch**, and follow the instructions.

#. Start developing as usual and keep in mind that regular commits
   make it easier to follow the development. Remember to update the
   author and date information in the files you modify.

#. Push changesets to the new branch::

    bzr push lp:<path-to-branch-location>

   The first time you push to this location you should use the
   ``--use-existing-dir`` option.

#. Once you have completed your work, you should propose it for
   merging into the DOLFIN main branch (via the Launchpad system). A
   developer with write access to the main branch will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for using branches for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

.. _contributing_patches:

Creating a patch
================

As an alternative to creating a branch, you may submit a patch. The
following instructions show how to create and submit a patch for
DOLFIN.

#. Get the current development branch::

    bzr branch lp:dolfin

#. Modify the files.

#. If your contribution consists of new files, add those to the
   branch::

    bzr add <files>

   where ``<files>`` is the list of new files. Do not add temporary or
   binary files. No action is necessary for previously existing files
   which have been modified.

#. Remember to update the author and date information in the files you
   modify.

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
   developer with write access to the main branch will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for creating a patch for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

.. _bzr_workflow:

Recommended Bazaar workflow
===========================

When working together with others on a code, it often happens that one
needs to merge changes from two or more branches of the same code. The
following is a recommended Bazaar workflow for handling merges. It
applies mainly to members of core teams that have write access to the
main branches, but may also be of use to others.

#. Create a shared repository for branches::

    bzr init-repo foo

   This creates a directory named ``foo`` which can hold several
   branches that share data, which not only saves disk space but also
   speeds up merging and branching.

#. Enter the shared repository::

    cd foo

#. Checkout the main branch of the project from Launchpad::

    bzr checkout lp:foo trunk

   This creates a *bound* branch of the project in the directory
   ``trunk``. Commits in this directory will result in a commit
   in the main Launchpad branch.

#. Create a branch for local work::

    bzr branch trunk work

#. Make any changes, commits, merges etc. inside the ``work``
   directory::

    cd work
    <work>
    <work>
    <work>
    bzr commit

#. When you want to transfer your changes to the main branch, first
   try to push your changes directly to the main branch::

    bzr push lp:foo

#. If that fails, which can happen if someone else has pushed changes
   to the main branch before you, a merge is necessary. The point now
   is that this merge should be carried out *from* the main
   branch. The merge should not be carried out inside the ``work``
   directory and then pushed to the main branch (as that will create a
   warning about revisions being removed from the main branch). Here's
   how to carry out the merge::

    cd ../trunk
    bzr update
    bzr merge ../work
    bzr commit -m "merge work on <stuff>"

   This will merge the changes made in ``work`` and transfer those
   changes to the main Launchpad branch.

   Some FEniCS projects have explicitly set the Bazaar flag
   ``append_revisions_only``, which will issue an error message if an
   attempt is made to push a merge from ``work``.
