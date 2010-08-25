.. Developer information.

.. _developer:


#####################
Developer information
#####################

This page contains information for FEniCS developers.  FEniCS development
tools are hosted Launchpad. This page explains how to obtain source
code from `Launchpad <https://launchpad.net/>`_, subscribe to developer
mailing lists and take part in the development of FEniCS. An overview
of all FEniCS projects on Launchpad can be found
`here <https://launchpad.net/fenics-project>`_.

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


