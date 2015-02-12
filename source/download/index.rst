.. How to download and install FEniCS projects.

.. _download:

.. include:: icons.rst

########
Download
########

.. _binary_packages:

*********************************
FEniCS versions and release notes
*********************************

The latest stable release of **FEniCS version 1.5.0** released
2015-01-12. For more information about FEniCS releases and changes,
take a look at :ref:`release_notes`. Note that binary packages will
appear a few days after the release of a new version.

***************
Binary packages
***************

This is the main FEniCS download, suitable for most users. It includes
everything needed to get you started with FEniCS.

.. raw:: html
    :file: index.inc

For information on user contributed binary packages, see
:ref:`contributed_packages`.

**********************************************
Installation from source (`fenics-install.sh`)
**********************************************

You may also choose to install FEniCS directly from source. This may be
done by running the following command:

    curl -s http://fenicsproject.org/fenics-install.sh | bash

Running this command will build a local installation of FEniCS, including
essential dependencies such as PETSc. The installation relies on
`HashDist <http://hashdist.github.io/>`__. Before running the script,
you may wish to download and expect its contents.

Requirements:

* A standard Unix environment (Linux or Mac)
* The Git version control system, available as the package ``git`` in most Linux and Mac package managers
* Mac users will also need to install Xcode from the Mac App Store

*************************
Nightly snapshot releases
*************************

Every night, FEniCS snapshot releases are automatically generated for
Ubuntu and Mac OS X. They are made available at our :ref:`snapshots page
<snapshot_releases>`.

*******************
Development version
*******************

For information on accessing the development repositories, see
:ref:`developers`.

***************
Data and meshes
***************

A collection of meshes for free use with FEniCS is available
:ref:`here <data>`.

***************
Troubleshooting
***************

Visit the :ref:`troubleshooting page <troubleshooting>` if you have
problems installing FEniCS.

.. toctree::
   :hidden:
   :glob:

   *
