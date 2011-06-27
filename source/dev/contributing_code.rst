.. _developers_contributing_code:

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

.. _before_committing:

Before committing your work
===========================

Before committing any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying. If you are
submitting code for a new feature, it is expected that the new feature
is accompanied by a suitable set of unit tests.

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
   :ref:`license <about_license>` section.

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

.. _bzr_branch_workflow:

Recommended Bazaar workflow
===========================

.. todo::

   Copy mailinglist discussion here.
