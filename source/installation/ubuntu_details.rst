########
Download
########

.. image:: images/ubuntu.png
    :align: left

`Click to install FEniCS <apt://fenics>`__

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
<https://launchpad.net/~fenics/+archive/ppa>`__ (PPA) to your Ubuntu
system. You can do this by running the following commands in a terminal:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics/ppa
    sudo apt-get update
    sudo apt-get install fenics

If you exeperience problems with the FEniCS packages after adding the
FEniCS PPA, make sure that all packages are up-to-date by running the
following command in a terminal::

    sudo apt-get dist-upgrade
