.. _ubuntu_details:

####################################
Installation instructions for Ubuntu
####################################

.. image:: images/debian.png
    :align: left

`Install FEniCS <apt://fenics>`__

FEniCS is included as part of Ubuntu GNU/Linux. To install, simply click
on the link above or run the following command in a terminal::

    sudo apt-get install fenics

.. warning::

    Since Ubuntu has a policy of not updating the version of a shipped
    application, the FEniCS version available in the standard Ubuntu
    archive might be out-of-date. To keep up with the latest FEniCS
    releases, add the Ubuntu PPA as described below.

.. _ubuntu_ppa:

**********
Ubuntu PPA
**********

To make sure you always have the latest FEniCS release, add the `FEniCS
Personal Package Archive
<https://launchpad.net/~fenics-packages/+archive/fenics>`__ (PPA) to
your Ubuntu system. You can do this by running the following commands in
a terminal:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics-packages/fenics
    sudo apt-get update
    sudo apt-get install fenics
    sudo apt-get dist-upgrade

.. note::

    The PPA will be deactivated if you later upgrade to a newer Ubuntu
    release, so you should run these commands again after a release
    upgrade.

================
Removing the PPA
================

If you want to remove the FEniCS PPA and restore everything back to the
default packages from the standard Ubuntu archive, run the following
command in a terminal:

.. code-block:: sh

    sudo ppa-purge ppa:fenics-packages/fenics

##############
Older releases
##############

Older versions of the Ubuntu packages are also :ref:`available
<older_releases>`.
