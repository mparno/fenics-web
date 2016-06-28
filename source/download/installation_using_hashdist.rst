.. _installation_using_hashdist:

***************************************
Installation from source using HashDist
***************************************


Requirements
============

Requirements for Linux (Ubuntu package name in parentheses):

* The Python 2.7 development files (``python-dev``)
* The Git version control system (``git``)
* A C compiler and a C++ compiler with C++11 support
  (``build-essential``)
* A Fortran compiler (``gfortran``)
* The OpenGL development files (``freeglut3-dev``)

Requirements for OS X:

* The latest version of Xcode `from the Apple developer website
  <https://developer.apple.com/downloads/index.action>`__ or get it
  `using the Mac App Store
  <http://itunes.apple.com/us/app/xcode/id497799835>`__
* The Xcode Command Line Developer Tools (run ``xcode-select --install`` in
  a terminal after installing Xcode)
* The Git version control system


Configuration
=============

The ``fenics-install.sh`` script is a nice tool to get a FEniCS
installation up and running with little effort. However, the default
installation given by ``fenics-install.sh`` is not always what we
want. Maybe we want to build without VTK or maybe we want to build our
own branch of DOLFIN. Maybe we need to set the paths to a specific
compiler or maybe we are on a cluster and want to use modules from the
modules system. The ``fenics-install.sh`` script can also be used in
these cases, but we have to modify the HashDist profile:

.. code-block:: bash

    git clone https://bitbucket.org/fenics-project/fenics-developer-tools.git
    cd fenics-developer-tools
    cp install/profiles/fenics.Linux.yaml fenics.yaml
    editor fenics.yaml  # modify the profile for your needs
    ./install/fenics-install.sh fenics.yaml

.. note::

    Use one of the other profiles under ``install/profiles`` as your
    template if you are not on Linux.

What the ``fenics-install.sh`` script actually does is to first
download HashDist and HashStack, then create a HashDist profile for
your platform and finally build the FEniCS stack using this
profile. Below, we will see how to do these steps manually and how to
modify the HashDist profile for our needs.


Installing HashDist
===================

`HashDist <https://hashdist.github.io>`__ can be installed using
``pip``::

    pip install --upgrade https://github.com/hashdist/hashdist/archive/master.zip

.. note::

    You might need to use ``sudo`` when running ``pip``.

HashDist can now be initialize by running::

    hit init-home

This will create the ``~/.hashdist`` directory and you might want to
modify ``~/.hashdist/config.yaml`` before you continue.


Creating the HashDist profile
=============================

First, clone a copy of `HashStack
<https://github.com/hashdist/hashstack>`__ from GitHub::

    git clone https://github.com/hashdist/hashstack.git

Enter the ``hashstack`` directory and create a new profile for
instance by copying one of the profiles from the ``examples``
directory::

    cd hashstack
    cp examples/fenics-1.6.0.linux.yaml fenics.yaml

We will in the following assume that the profile is called
``fenics.yaml``.

Now you can start modifying the profile for your needs. The `HashDist
documentation
<http://hashdist.readthedocs.org/en/latest/specs.html>`__ contains
information about how to create a profile. Here we will only include
some examples that could be useful when building FEniCS.


###########################
Adding or removing packages
###########################

Adding or removing a package is easy. Simply add or remove it from the
``packages`` section in your profile. Let us say you want to build
``libadjoint`` and ``dolfin-adjoint``. Simply add them to the
``packages`` list::

    packages:
      ...
      libadjoint:
        build_with: python
      dolfin-adjoint:

If a package you need is not avilable in HashStack, please consider
adding it yourself and make a `pull request
<https://github.com/hashdist/hashstack/compare>`__, or suggest that it
is added by `registering an issue
<https://github.com/hashdist/hashstack/issues/new>`__.


#########################
Using alternative sources
#########################

Sometimes we might want to use an alternative source for a
package. For example, assume that we want to use the latest
development version of DOLFIN. One way to do this is by running::

    hit fetchgit -p dolfin https://bitbucket.org/fenics-project/dolfin.git master

The output will be the latest git commit id (something like
``git:19b84ef42e9926450e8b3751adffd6b47078f322``), which can then be
used to override the source for DOLFIN in the profile::

    packages:
      ...
      dolfin:
        sources:
        - key: git:19b84ef42e9926450e8b3751adffd6b47078f322
	  url: https://bitbucket.org/fenics-project/dolfin.git

Similarly, if we want to build our own fork of DOLFIN, we can do::

    hit fetchgit -p dolfin https://bitbucket.org/<username>/dolfin.git master

The output can then be used to override the source as shown above.


#############################
Setting environment variables
#############################

HashDist overrides the ``PATH`` environment variable. If you for
instance want to use a compiler that is not available in the standard
path, then you need to add this to the ``parameters`` section in your
profile. Let us assume that we want to use a compiler that is
installed in ``/opt/bin``. This can be done by adding the following to
our profile::

    parameters:
      PATH: /opt/bin:/usr/bin:/bin

Similarily, if we are on a system that uses the Environment Modules
system, as used on many clusters, we can add commands such as ``module
load gcc/4.9.2`` to the ``PROLOGUE`` parameter. The commands in the
``PROLOGUE`` parameter are executed before building every package. You
can find examples of this in the HashStack repository (see
``fenics.abel.gnu.yaml`` and ``fenics.scinet.gnu.yaml``).

.. note::

    Variables defined in the ``parameters`` section should also be
    defined when you run your FEniCS programs.


Building the profile
====================

When you are satisfied with your profile, run the following command to
start building::

    hit build -j4 fenics.yaml

The ``-j`` option is used to build in multiple threads to speed up the
build. Here we build with 4 threads, but this can be adjusted as you
like.

When the build has finished, you will have a new directory ``fenics``
(the profile name without ``.yaml``), which contains the complete
installation of FEniCS. Before you can start running your FEniCS
programs, you need to set some environment variables:

.. code-block:: sh

    export PATH=$HOME/hashstack/fenics/bin:$PATH
    export CMAKE_PREFIX_PATH=$HOME/hashstack/fenics

Here we assume that the ``hashstack`` directory is located directly
under ``$HOME``.
