.. Developer information.

.. _developer:


#####################
Developer information
#####################

FEniCS development takes place on Launchpad. This page explains how to obtain
FEniCS source code from Launchpad, subscribe to developer mailing lists and
take part in the development of FEniCS.

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

bzr branch lp:project-name

For example, to obtain the source code for DOLFIN, run the command

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


