.. How to download and install FEniCS projects.

.. _download:

.. |dockerimage| image:: images/docker.png
.. |debianimage| image:: images/debian.png
.. |ubuntuimage| image:: images/ubuntu.png
.. |toolsimage|  image:: images/tools.png
.. |hspace| raw:: html

  &nbsp;

################################
FEniCS Download and Installation
################################

The easiest and best way to install FEniCS is to use our prebuilt
high-performance `Docker <https://www.docker.com>`_ images. FEniCS can
also be installed via package managers or compiled directly from
source. Follow the instructions below to get started with FEniCS.

*****************************************************************
|dockerimage| |hspace| Docker images (all platforms and versions)
*****************************************************************

`Docker <https://www.docker.com>`_ allows us to build and ship
consistent high-performance FEniCS installations for almost any
platform. To get started, follow these 3 simple steps:

Step 1: Install Docker itself. Mac and Windows users should install the
`Docker Toolbox <https://www.docker.com/products/docker-toolbox>`_
(this is a simple one-click install) and Linux users should
`follow these simple instructions
<https://docs.docker.com/linux/step_one/>`_.

Step 2: If running on Mac or Windows, make sure you have the
Docker Quickstart Terminal running. (This may take a little while
the first time so just be patient.)

Step 3: Run the prebuilt FEniCS image by the following simple command::

    docker run -ti quay.io/fenicsproject/stable

The first time you run this command, Docker will automatically fetch
the image for the latest stable version of FEniCS. For the latest
development version of FEniCS, just change ``stable`` to ``dev``.

To share files between your operating system and the FEniCS Docker
image, simply add the ``-v`` flag to tell Docker where your files
are, for example::

    docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable

For more details and tips on how to work with FEniCS and Docker, see
our `FEniCS Docker page
<http://fenics-containers.readthedocs.org/en/latest/>`_.

*******************************************************
|ubuntuimage| |hspace| Ubuntu packages (stable release)
*******************************************************

FEniCS is part of the Ubuntu (and Debian) GNU/Linux operating systems.
To install FEniCS, run the following command::

    sudo apt-get install fenics

For more details and tips on how to work with FEniCS in Ubuntu, see
our :ref:`FEniCS Ubuntu page <ubuntu_details>`.

*************************************************************************
|toolsimage| |hspace| Manual installation from source (all platforms and versions)
*************************************************************************

FEniCS can be built manually from source using standard installation
mechanisms for Python
(`Setuptools <https://pypi.python.org/pypi/setuptools>`_)
and C++
(`CMake <https://cmake.org/>`_).

The
`FEniCS source code
<https://bitbucket.org/account/user/fenics-project/projects/CORE>`_
is hosted in
`Git <https://git-scm.com/>`_ repositories on
`Bitbucket <https://bitbucket.org/>`_.

****************************************************************************
|toolsimage| |hspace| Automatic installation from source (all platforms and versions)
****************************************************************************

FEniCS can be built automatically from source via
`HashDist <http://hashdist.github.io/>`__. To build FEniCS, run the
following command::

    curl -s http://fenicsproject.org/fenics-install.sh | bash

Running this command will build a local installation of FEniCS.
Before running the script, you may wish to download and inspect its
contents.

For more details and tips on how to work with FEniCS in HashDist, see
our :ref:`FEniCS HashDist page <installation_using_hashdist>`.

********************
Contributed packages
********************

FEniCS is also available through a number of alternative package
managers.

For information on user contributed binary packages, see
our :ref:`FEniCS contributed packages page <contributed_packages>`.

****************************
FEniCS versions and releases
****************************

The latest stable release of FEniCS is version
1.6.0,
which was released on
2015-07-28.

For more information about FEniCS releases and
changes, take a look at :ref:`release_notes`. For information on
accessing the development repositories, see :ref:`developers`.

***************
Data and meshes
***************

A collection of meshes for free use with FEniCS is available
:ref:`here <data>`.

.. toctree::
   :hidden:
   :glob:

   *
