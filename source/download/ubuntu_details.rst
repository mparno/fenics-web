.. _ubuntu_details:

########
Download
########

.. image:: images/ubuntu.png
    :align: left

`Install FEniCS <apt://fenics>`__

#########################
Installation instructions
#########################

FEniCS is included as part of Ubuntu GNU/Linux (starting with
10.04/Lucid). To install, simply click on the link above or install from
the Ubuntu Software Center. FEniCS can also be installed by running the
following command in a terminal::

    sudo apt-get install fenics

**********
Ubuntu PPA
**********

Since Ubuntu has a policy of not updating the version of a shipped
application, the FEniCS version available in the standard Ubuntu archive
might be old. To keep up with the latest FEniCS releases, add the
`FEniCS Personal Package Archive
<https://launchpad.net/~fenics-packages/+archive/fenics>`__ (PPA) to
your Ubuntu system. You can do this by running the following commands in
a terminal:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics-packages/fenics
    sudo apt-get update
    sudo apt-get install fenics
    sudo apt-get dist-upgrade

.. warning::

    Ubuntu 10.04 LTS and 10.10 users should note that when adding the
    FEniCS PPA and installing the FEniCS packages, the Boost packages on
    the system will be upgraded to version 1.42. In most cases this is
    fine since only the Boost `development` packages (as in
    ``libboost-foo-dev``) will be replaced. However, if you want to keep
    the default Boost packages, then you should not add this PPA.

================
Removing the PPA
================

If you want to remove the FEniCS PPA and restore everything back to the
default packages from the standard Ubuntu archive, run the following
command in a terminal:

.. code-block:: sh

    sudo ppa-purge ppa:fenics-packages/fenics
