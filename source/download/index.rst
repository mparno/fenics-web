.. How to download and install FEniCS projects.

.. _download:

.. |dockerimage| image:: images/docker.png
.. |debianimage| image:: images/debian.png
.. |ubuntuimage| image:: images/ubuntu.png
.. |toolsimage|  image:: images/tools.png
.. |hspace| raw:: html

  &nbsp;

################################
FEniCS download and installation
################################

The easiest way to install FEniCS is to use our prebuilt,
high-performance `Docker <https://www.docker.com>`_ images. FEniCS can
also be installed via package managers or compiled directly from
source. Follow the instructions on this page to get started with FEniCS.

The latest stable release of FEniCS is version **2016.2.0**, which was
released on 2016-11-30. For more information about FEniCS releases and
changes, take a look at :ref:`release_notes`. For information on
accessing the development repositories, see :ref:`developers`.

*****************************************************************
Docker images (all platforms and versions) |hspace| |dockerimage|
*****************************************************************

`Docker <https://www.docker.com>`_ allows us to build and ship
consistent high-performance FEniCS installations for almost any
platform. To get started, follow these 2 steps:

#. Install Docker for your platform:
   `Windows <https://www.docker.com/products/docker-toolbox>`_,
   `Mac <https://www.docker.com/products/docker#/mac>`_ or
   `Linux <https://www.docker.com/products/docker#/linux>`_. 
#. Install the FEniCS Docker script::

    curl -s https://get.fenicsproject.org | bash

Notes:

* Windows users must (for now) continue to use the Docker Toolbox (*not* Docker
  for Windows) if they want to use the ``fenicsproject`` script described below.
* If using the Docker Toolbox (macOS versions < 10.10 or Windows), make sure
  you run all commands inside the Docker Quickstart Terminal.

Once both Docker and the FEniCS Docker script have been installed, you
can easily start a FEniCS session by running the following command::

    fenicsproject run

The FEniCS Docker script can also be used to create persistent sessions
(``fenicsproject create myproject`` followed by ``fenicsproject run
myproject``) or to run different versions of FEniCS (``fenicsproject run
dev``). To see more options, run the following command::

    fenicsproject help

For all ``fenicsproject`` commands, the contents of the current working
directory will be shared into the project at ``~/shared``.

We regularly release new Docker images with updated versions of FEniCS
and its supporting libraries. To upgrade to the latest version, run::

    fenicsproject pull

To upgrade to the latest version of the ``fenicsproject`` script, run::

    curl -s https://get.fenicsproject.org | bash

For more details and tips on how to work with FEniCS and Docker, see
our `FEniCS Docker page
<http://fenics-containers.readthedocs.org/en/latest/>`_.

*******************************************************
Ubuntu packages (stable release) |hspace| |ubuntuimage|
*******************************************************

To install FEniCS on Ubuntu, run the following commands::

    sudo add-apt-repository ppa:fenics-packages/fenics
    sudo apt-get update
    sudo apt-get install fenics
    sudo apt-get dist-upgrade

This will add our `PPA <https://launchpad.net/ubuntu/+ppas>`_ for
FEniCS to your package sources and install the latest stable version
of FEniCS. Note that FEniCS is also available from the official
Ubuntu (and Debian) repositories but may outdated, depending on
which release of Ubuntu you are running.

For more details and tips on how to work with FEniCS in Ubuntu, see
our :ref:`FEniCS Ubuntu page <ubuntu_details>`.


**********************************************************************************
Manual installation from source (all platforms and versions) |toolsimage| |hspace|
**********************************************************************************

FEniCS can be built manually from source using standard installation
mechanisms for Python (`Setuptools
<https://pypi.python.org/pypi/setuptools>`_) and C++ (`CMake
<https://cmake.org/>`_).

The `FEniCS source code
<https://bitbucket.org/account/user/fenics-project/projects/CORE>`_ is
hosted in `Git <https://git-scm.com/>`_ repositories on `Bitbucket
<https://bitbucket.org/>`_.


*************************************************************************************
Automatic installation from source (all platforms and versions) |toolsimage| |hspace|
*************************************************************************************

FEniCS can be built automatically from source via
`HashDist <http://hashdist.github.io/>`__. To build FEniCS, run the
following command::

    curl -s https://fenicsproject.org/fenics-install.sh | bash

Running this command will build a local installation of FEniCS.
Before running the script, you may wish to download and inspect its
contents.

For more details and tips on how to work with FEniCS in HashDist, see
our :ref:`FEniCS HashDist page <installation_using_hashdist>`.


********************
Contributed packages
********************

FEniCS is also available through a number of alternative package
managers. For information on user contributed binary packages, see our
:ref:`FEniCS contributed packages page <contributed_packages>`.


***************
Data and meshes
***************

A collection of meshes for free use with FEniCS is available
:ref:`here <data>`.

.. toctree::
   :hidden:
   :glob:

   *
