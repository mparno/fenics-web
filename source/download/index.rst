.. How to download and install FEniCS projects.

.. _download:

.. include:: icons.rst

#################################
FEniCS versions and release notes
#################################

FEniCS can be installed from binary files or compiled from source
code. Different installation methods are described on this page.
Visit the :ref:`troubleshooting page <troubleshooting>` if you have
problems installing FEniCS.

The latest stable release of FEniCS is version 1.6.0, which was on
released 2015-07-28. For more information about FEniCS releases and
changes, take a look at :ref:`release_notes`.  For information on
accessing the development repositories, see :ref:`developers`.


***************
Binary packages
***************

This is the main FEniCS download, suitable for most users. It includes
everything needed to get you started with FEniCS.

.. raw:: html
    :file: index.inc

For information on user contributed binary packages, see
:ref:`contributed_packages`.

*************
Docker images
*************

Docker_ allows us to build, ship and run FEniCS on almost any platform with just
one command::

    docker run -ti fenicsproject/stable

You will be presented with an environment where you can run the latest version
of FEniCS with all of the complex dependencies (HDF5, PETSc, SLEPc, petsc4py,
slepc4py, MPI etc.) already compiled.

To get started, install Docker_ for your platform (Windows, Mac OS X, Linux,
cloud environment, and many others) using the instructions at
https://docs.docker.com/engine/installation/. There are no other installation
requirements or dependencies.

We provide a range of Docker images to meet the varied needs of developers,
advanced and first time users. For more details see the instructions at
https://bitbucket.org/fenics-project/docker.

.. _Docker: https://docker.com

**********************************
Installation from *conda* packages
**********************************

Whether you are already using the Anaconda_ Python distribution from Continuum
Analytics or just their conda_ package manager in Linux you can already install
FEniCS precompiled binaries. They are not part of the official Anaconda
installation but they can be downloaded from a binstar_ channel::

    conda create --name fenics27 python=2.7
    source activate fenics27
    conda install fenics --channel juanlu001

For troubleshooting and alternative options for older systems, see
:ref:`installation_using_conda`.

.. note::

    These packages are provided by Juan Luis Cano and their sources
    can be found at https://github.com/juanlu001/fenics-recipes.

.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _conda: http://conda.io/
.. _binstar: https://binstar.org/


**********************
Virtual machine images
**********************

A virtual machine image that includes the most recent FEniCS release
is available at
http://fenicsproject.org/pub/virtual/fenics-latest.ova. The image can
be run using a virtual machine manager, such as `VirtualBox
<https://www.virtualbox.org/>`_. The username and password for the
virtual machine are both ``fenics``.

The virtual machine image is recommended for systems for which a
binary installer is not available. The image is particularly suitable
when a consistent FEniCS environment across systems is required, such
as courses using FEniCS.

The image requires a 64-bit host operating system.


**********************************************
Installation from source (`fenics-install.sh`)
**********************************************

You may also choose to install FEniCS directly from source. This may be
done by running the following command::

    curl -s http://fenicsproject.org/fenics-install.sh | bash

Running this command will build a local installation of FEniCS, including
essential dependencies such as PETSc. The installation relies on
`HashDist <http://hashdist.github.io/>`__. Before running the script,
you may wish to download and inspect its contents.

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
  
See :ref:`installation_using_hashdist` if you need to tweak the build.


***************
Data and meshes
***************

A collection of meshes for free use with FEniCS is available
:ref:`here <data>`.



.. toctree::
   :hidden:
   :glob:

   *
