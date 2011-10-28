.. _ubuntu_details:

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

Ubuntu 10.04 LTS and Ubuntu 10.10 users will also need to add the
`FEniCS experimental builds PPA
<https://launchpad.net/~fenics-packages/+archive/fenics-exp>`__ to get
the latest FEniCS releases. Use the following set of commands instead of
the above:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics/ppa
    sudo add-apt-repository ppa:fenics-packages/fenics-exp
    sudo apt-get update
    sudo apt-get install fenics

.. warning::

    Adding the `FEniCS experimental builds PPA
    <https://launchpad.net/~fenics-packages/+archive/fenics-exp>`__ will
    also upgrade the Boost packages on your system to version 1.42. In
    most cases this is unproblematic since only the Boost `development`
    packages (as in ``libboost-foo-dev``) will be replaced. However, if
    you want to keep the default Boost packages, then do not add this
    PPA.

    If you accidentally added this PPA and want to go back to the
    default Boost packages, running ``sudo ppa-purge
    ppa:fenics-packages/fenics-exp`` should do the trick.

.. note::

    If you experience problems with the FEniCS packages after adding the
    FEniCS PPA, make sure that all packages are up-to-date by running
    ``sudo apt-get dist-upgrade`` in a terminal.
