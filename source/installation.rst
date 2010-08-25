.. How to download and install FEniCS projects.

.. _installation:

############
Installation
############

Obtaining and installing FEniCS has never been easier. Packages are
available for Debian and Ubuntu users and we provide a framework
(`Dorsal <http://www.fenics.org/dorsal/>`_) for simple builds of
FEniCS on a multitude of platforms.

**********************
Debian/Ubuntu packages
**********************

Debian
======

FEniCS is included as part of Debian GNU/Linux
(testing/squeeze and unstable/sid). To install, simply click `here <apt://fenics>`_
(requires ``aptlinex``). FEniCS can also be installed by running the
following command::

    sudo apt-get install fenics

Ubuntu
======

FEniCS is included as part of Ubuntu GNU/Linux (starting with
10.04/Lucid). To install, simply click `here <apt://fenics>`_
(requires ``apturl``) or install from the Ubuntu Software
Center. FEniCS can also be installed by running the following
command::

    sudo apt-get install fenics

Ubuntu PPA
----------

The official Ubuntu repositories may not always contain the latest
stable FEniCS releases. To stay current with the latest releases, add
the `FEniCS Personal Package Archive
<https://launchpad.net/~fenics/+archive/ppa>`_ (PPA) to your Ubuntu
system. You can do this by going to **System > Administration >
Software Sources > Third-Party Software** and clicking on
**Add**. Then type in ``ppa:fenics/ppa`` and click on **Add Source**
and then **Close**. A dialog box should appear, asking whether you'd
like to update the list of repositories. Select **Reload** to update
the list. You can then install FEniCS from the Ubuntu Software Center.

Here is a list of commands for installing FEniCS from the FEniCS PPA
for those preferring the command-line:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics/ppa
    sudo apt-get update
    sudo apt-get install fenics

.. note::

    The ``add-apt-repository`` command is not available on older (pre
    9.10) Ubuntu systems. Please see the `FEniCS PPA page
    <https://launchpad.net/~fenics/+archive/ppa>`_ for instructions on
    these systems.

***********************************
Automated installation using Dorsal
***********************************

If you are not running a Debian or Ubuntu operating system, you need
to download and compile FEniCS and all its dependencies. Luckily, this
can be handled easily using `Dorsal`_.
Dorsal is a simple shell script that automates the process of
fetching, compiling, and installing the various FEniCS sub-projects
and their requisite dependencies on Unix-like operating systems.

The following platforms are officially supported by Dorsal:

* Debian GNU/Linux testing, unstable
* Fedora 13
* Gentoo (and Sabayon) Linux
* Mac OS X (Snow Leopard)
* openSUSE 11.3
* Red Hat Enterprise Linux 6
* Ubuntu 9.10, 10.04 LTS, 10.10
* Ubuntu (10.04/Lucid and 10.10/Maverick)

A number of other platforms are supported through user contributions.
It is easy to extend support to other similar platforms. Please let us
know if you are interested, and we can work together toward adding
this functionality.

In order to install FEniCS using Dorsal, simply follow these steps:

#. Fetch the most recent version of Dorsal from its
   `download page <https://launchpad.net/dorsal/+download>`_.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and tweak dorsal.cfg to your liking.
#. Invoke the build script by running::

    ./dorsal.sh

#. At this point, Dorsal attempts to guess your operating system
   (platform) and provides a list of instructions to ensure that your
   system is ready for installing FEniCS. Go through these
   instructions carefully, and copy and paste appropriate commands
   into a separate terminal window to prepare your system.
#. Once you have completed these steps, hit enter to begin the
   installation!

Once the build starts, you should see a lot of compilation-related
text scrolling past your screen and eventually end up with a complete,
up-to-date installation of FEniCS. It will take some time to build all
the libraries, so be patient.

Don't forget to follow any post-build instructions before rushing off
to try the demos!

*******************************
Manual installation from source
*******************************

You can also build and install FEniCS components manually from the source code.
You will need to install the Python packages
`FFC <http://www.fenics.org/ffc>`_,
`FIAT <http://www.fenics.org/fiat>`_,
`Instant <http://www.fenics.org/instant>`_,
`Viper <http://www.fenics.org/fenics-viper>`_, and
`UFL <http://www.fenics.org/ufl>`_.
You will also need to install the C++/Python packages
`UFC <http://www.fenics.org/ufc>`_ and
`DOLFIN <http://www.fenics.org/dolfin>`_.
Instructions for how to install these packages are given below.

Python packages (FFC, FIAT, Instant, Viper and UFL)
===================================================

#. First, `download the source code <http://www.fenics.org/wiki/Download#Download_the_source_code>`_
   for FFC, FIAT, Instant, Viper and UFL. Then for each of these packages,
   follow the instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following command::

    sudo python setup.py install

This will install the packages system wide on your computer. You may
easily change the installation directory. For example, you may wish to
install the packages in a directory named ``local`` in you home
directory. This may be accomplished by running the following command::

    sudo python setup.py install --prefix=~/local

C++/Python packages (UFC and DOLFIN)
====================================================

#. First, `download the source code`_
   for UFC and DOLFIN. Then for each of these packages, follow the
   instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following commands::

    cmake .
    make
    sudo make install

This will install the packages system wide on your computer. You may
easily change the installation directory. For example, you may wish to
install the packages in a directory named ``local`` in you home
directory. This may be accomplished by replacing the first of the
above commands by the following command::

    cmake . -DCMAKE_INSTALL_PREFIX=~/local

You may also want to consider using a graphical frontend for CMake
like either ``cmake-gui`` or ``ccmake`` which both provide a simple
way to configure the installation.

During the configuration phase of DOLFIN (calling ``cmake``,
``cmake-gui``, or ``ccmake``), you will be notified of any missing
dependencies. If a required package is missing, you will need to
install that package and configure DOLFIN again. If an optional
package is missing, you may choose to continue with the installation
but some functionality may be missing.

***************
Troubleshooting
***************

If for some reason you were unable to obtain, build, or install a
FEniCS component, the best place to start is the `Launchpad Answers
<https://help.launchpad.net/Answers>`_ page for the project which
causes the problem:

* `DOLFIN <http://answers.launchpad.net/dolfin>__
* `Dorsal <https://answers.launchpad.net/dorsal>`__
* `FFC <https://answers.launchpad.net/ffc>`__
* `FIAT <https://answers.launchpad.net/fiat>`__
* `Instant <https://answers.launchpad.net/instant>`__
* `Viper <https://answers.launchpad.net/fenics-viper>`__
* `UFC <https://answers.launchpad.net/ufc>`__
* `UFL <https://answers.launchpad.net/ufl>`__

If all else fails, you can always try the fenics@lists.launchpad.net
mailing list.

If you have suggestions for improving the information found on this
page with respect to the installation procedure of FEniCS, you are
more than welcome to file a bug report or register a blueprint on
`FEniCS Documentation <https://launchpad.net/fenics-doc>`_.

