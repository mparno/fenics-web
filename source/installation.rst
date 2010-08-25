.. How to download and install FEniCS projects.

.. _installation:

############
Installation
############

Obtaining and installing FEniCS software has never been easier.
Packages are available for Debian and Ubuntu users and Dorsal provides a
framework for custom builds of FEniCS software and dependencies.
If for some reason you were unable to obtain, build or install a FEniCS
component, the best place to start is the
`Launchpad Answers <https://help.launchpad.net/Answers>`_ page for the project
which cause the problem:

* `DOLFIN <https://answers.launchpad.net/dolfin>`_
* `Dorsal <https://answers.launchpad.net/dorsal>`__
* `FFC <https://answers.launchpad.net/ffc>`_
* `UFC <https://answers.launchpad.net/ufc>`_
* `UFL <https://answers.launchpad.net/ufl>`_

(Visit the `FEniCS Project <https://launchpad.net/fenics-project>`_ page for a
complete list of FEniCS components)

If everything else fails you can always try the fenics@lists.launchpad.net
mailing list.
If you have suggestions for improving the information found on this page with
respect to the installation procedure of FEniCS, you are more than welcome to
file a bug report or register a blueprint on
`FEniCS Documentation <https://launchpad.net/fenics-doc>`_.

**********************
Debian/Ubuntu packages
**********************

Debian
======

The latest stable FEniCS releases are available in Debian sid. To
install, simply click `here <apt://fenics>`_ (requires ``aptlinex``) or
run the following command as root::

    aptitude install fenics

Ubuntu
======

FEniCS is available in the official Ubuntu repositories as of Ubuntu
10.04 (Lucid). Simply click `here <apt://fenics>`_ to install
(requires ``apturl``) or install from the Ubuntu Software
Center. FEniCS can also be installed by running the following
command::

    sudo apt-get install fenics

FEniCS PPA
--------------

The official Ubuntu repositories may not always contain the latest
stable FEniCS releases. To stay current with the latest releases, add
the `FEniCS Personal Package Archive
<https://launchpad.net/~fenics/+archive/ppa>`_ (PPA) to your Ubuntu
system. You can do this by going to **System > Administration >
Software Sources > Third-Party Software** and click on **Add**. Then
type in ``ppa:fenics/ppa`` and click on **Add Source** and then
**Close**. A dialog box should appear, asking whether you'd like to
update the list of repositories. Select **Reload** to update the
list. You can then install FEniCS from the Software Center.

Here is a list of commands for installing FEniCS from the FEniCS PPA
for those preferring the command line:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics/ppa
    sudo apt-get update
    sudo apt-get install fenics

.. note::

    The ``add-apt-repository`` command is not available on older (pre
    9.10) Ubuntu systems. Please see the `FEniCS PPA page
    <https://launchpad.net/~fenics/+archive/ppa>`_ for instructions on
    these systems.

************
Using Dorsal
************

The easiest way to install FEniCS on a UNIX-like operating environment
is to use `Dorsal <https://launchpad.net/dorsal>`_.
Dorsal is a simple shell script that automates the process of fetching,
compiling and installing the various FEniCS sub-projects and their requisite
dependencies.
It currently supports officially the following platforms:

* Debian GNU/Linux 5.0, unstable
* Fedora 13
* Mac OS X 10.5, 10.6
* openSUSE 11.3
* Red Hat Enterprise Linux 5
* Ubuntu 9.10, 10.04 LTS, 10.10
* Gentoo (and Sabayon) Linux

A number of other platforms are supported through user contributions.
It is easy to extend support to other similar platforms. Please let us
know if that you are interested, and we can work together toward
adding this functionality.

In order to install the FEniCS project on one of the supported
platforms listed above, follow the following steps:

#. Fetch the most recent version of Dorsal from its
   `download page <https://launchpad.net/dorsal/+download>`_.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and tweak dorsal.cfg to your liking.
#. Invoke the build script by running::

    $ ./dorsal.sh

#. At this point, Dorsal attempts to guess your operating environment
   (platform) and provides a list of instructions to ensure your system is
   ready for installing FEniCS. Go through these carefully, and copy and paste
   appropriate commands into another terminal window to prepare your system.
#. Once you've completed these steps, hit enter to begin the installation!


If all goes according to plan, you should see a lot of
compilation-related text scroll past your screen and eventually end up
with a complete, up-to-date installation of various FEniCS projects
and their corresponding dependencies. It will take some time to build all the
libraries, so be patient.

Don't forget to follow any post-build instructions before rushing off to try
the demos!


*******************
Manual installation
*******************

If you wish to  install FEniCS components manually you can easily do so.

Python packages (FFC, FIAT, Viper and UFL)
==========================================

For a system wide installation of the development version of FFC do::

    $ bzr branch lp:ffc
    $ cd ffc
    $ sudo python setup.py install

Alternatively, you can specify the installation path, if for instance you don't
have root privileges, for a local installation do::

    $ python setup.py install --prefix=/home/user/local

Installing UFC and UFL follows the same procedure only the project name ``ffc``
must be ``ufc`` and ``ufl`` respectively.

To install previous stable releases visit:

* `download FFC <https://launchpad.net/ffc/+download>`_
* `download UFC <https://launchpad.net/ufc/+download>`_
* `download UFL <https://launchpad.net/ufl/+download>`_

download the desired tar ball, unpack and install using the same procedure as
outlined above.

UFC
===


DOLFIN
======

For a system wide installation of the development version of DOLFIN do::

    $ bzr branch lp:dolfin
    $ cd dolfin
    $ scons configure <options>

Run the command::

    $ scons configure --help

to see available ``<options>``.
Then build and install DOLFIN by running the commands::

    $ scons
    $ scons install

You can also install DOLFIN in the local DOLFIN tree simply by running the
script::

    $ ./scons.local

Visit `download DOLFIN <https://launchpad.net/dolfin/+download>`_ to get
previous stable releases; unpack and install following the above procedure.

