.. _developers_getting_code:

*************************
Obtaining the source code
*************************

FEniCS uses `Bitbucket <http://bitbucket.org>`__ for hosting code. Each FEniCS component has a
`Git <http://git-scm.com/>`_ repository on Bitbucket that contains all
source code (including the entire development history). The
repositories are readable for everyone, but write access is only
granted to the members of the core teams.

Accessing the development repositories
======================================

To access the development repositories, you first need to install the
revision control system Git. Visit the `Git web page
<http://git-scm.com/>`__ for instructions on how to install Git on
your operating system.

Once Git has been installed, you can access the development
repository of any FEniCS project by the ``git`` command. For example,
to check out the source code for DOLFIN, simply issue the following
command::

    git clone git@bitbucket.org:fenics-project/dolfin.git

Notifications of updates
========================

Developers should subscribe to notifications of changes made to the
source code by visiting the repository on Bitbucket and clicking the
'follow' button.

Links to the source repositories for FEniCS projects can be found on the
`FEniCS Bitbucket page <https://bitbucket.org/fenics-project>`__.

FEniCS Developer Tools
======================

Developers should take a look at the
`FEniCS Developer Tools <https://bitbucket.org/fenics-project/fenics-developer-tools>`__
repository. This contains scripts that are very useful for developers in setting
up a good development environment for FEniCS. In particular, consider using the scripts
``fenics-install-all.sh`` and ``fenics-install-component.sh``.
