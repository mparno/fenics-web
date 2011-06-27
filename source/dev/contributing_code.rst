.. _developers_contributing_code:

*****************
Contributing code
*****************

The main repository, or branch, for each FEniCS component is owned by
the core team of that component. Therefore, unless you are a member of
the core team, you will not be able to upload any code to the main
repository directly. Instead you will have to create your own
:ref:`branch <contributing_branches>` or submit a :ref:`patch
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

.. note::
    Need to write something here about needing to submit the license
    agreement.

.. _contributing_branches:

Creating a branch
=================

A simple way to submit your changes is to create a branch on Launchpad
and submit a merge request. This allows the maintainers to pull the
source code from your branch, review it and then push to the main
development repository. Using DOLFIN as an example, simply do the
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
   developer with write access to the main repository will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for using branches for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

.. _contributing_patches:

Creating a patch
================

As an alternative to creating a branch, you may submit a patch.  The
following instructions show how to create and submit a patch for
DOLFIN.

#. Get the current development branch::

    bzr branch lp:dolfin

#. Modify the files.

#. If your contribution consists of new files, add those to the
   repository::

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
   developer with write access to the main repository will then review
   the code and merge it with the main branch (assuming that it passes
   the code review).

The procedure for creating a patch for other FEniCS components is
identical (with ``dolfin`` replaced by the relevant component name).

.. _bzr_branch_workflow:

Recommended Bazaar workflow
===========================

.. note::
    Need to write about Bazaar usage here.
