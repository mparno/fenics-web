.. How to download and install FEniCS projects.

.. _installation:

############
Installation
############

Obtaining and installing FEniCS have never been easier. Packages are
available for Debian/Ubuntu, Mac OS X and Windows users, and we provide
a framework (`Dorsal <http://launchpad.net/dorsal/>`_) for
simple builds of FEniCS on a multitude of platforms. It is also possible
to download source tarballs for all components. For information on
accessing the development repositories, see :ref:`developer`.


*****************
Prebuilt binaries
*****************

Debian
======

FEniCS is included as part of Debian GNU/Linux (squeeze/testing and
sid/unstable). To install, simply click `here <apt://fenics>`_
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

Ubuntu has a policy of not updating the version of a shipped application,
so the FEniCS version available might be old. To keep up with the latest
releases, add the `FEniCS Personal Package Archive
<https://launchpad.net/~fenics/+archive/ppa>`_ (PPA) to your Ubuntu
system. You can do this by going to **System > Administration >
Software Sources > Other Software** and click on
**Add**. Then type in ``ppa:fenics/ppa`` and click on **Add Source**
and then **Close**. A dialog box should appear, asking whether you
would like to update the list of repositories. Select **Reload** to
update the list. You can then install FEniCS from the Ubuntu Software
Center.

Here is a list of commands for installing FEniCS from the FEniCS PPA
for those preferring the command-line:

.. code-block:: sh

    sudo add-apt-repository ppa:fenics/ppa
    sudo apt-get update
    sudo apt-get install fenics

.. note::

    Make sure to upgrade your packages after adding the FEniCS PPA. This
    can be done simply by running ``sudo apt-get dist-upgrade`` in a
    terminal or by going to **System > Administration > Update Manager**
    and upgrade the packages from there.

Mac OS X
========

The FEniCS binary for OS X runs on Intel 10.6 only. *10.5 and older and
the PowerPC architecture are not supported.* XCode 3 (available on the
OS X install disc) is required to run the binary.

* Install `FEniCS <http://www.fenicsproject.org/pub/software/fenics/fenics-11.05-osx10.6.dmg>`_

The installer should be mounted automatically after the download
completes. Otherwise, double-click on the ``.dmg`` file to fire up the
installer. Then simply drag the FEniCS icon into the Applications folder
to complete the installation.

To start using FEniCS, click on the FEniCS icon in the Applications
folder. This will bring up a terminal with everything set up to work
with FEniCS. Alternatively, source the file
``/Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf``
to set up the necessary paths to work with FEniCS. You can do this by
adding the following line to the ``.profile`` file in your home
directory::

    source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf

This will make FEniCS available whenever you start a new terminal.

.. note::

    If the FEniCS icon does not show up in the Applications folder
    after the installation, run ``killall -KILL Dock`` in a terminal to
    restart the Dock.

.. warning::

    This binary package is built against system Python and will not work
    with MacPorts Python or similar.

Windows
=======

The Windows binary installer will install everything needed to run
FEniCS on Windows, including MinGW compilers, Python, CMake, SWIG, and
others. It has been verified to run on XP, Vista, and on Windows 7.

You can download the installer by clicking `here
<http://www.fenicsproject.org/pub/software/fenics/fenics-11.05-mingw32.exe>`_.
Then double-click on the file and follow the instructions. When the
installation is complete, simply start the FEniCS command shell from the
start menu and you are ready to start working.

.. note::

    The Windows installer does not currently include PETSc or Trilinos
    which makes it unsuitable for solving any real problems.

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

* Debian GNU/Linux (squeeze/testing, sid/unstable)
* Fedora 13
* Gentoo Linux
* Mac OS X (Snow Leopard)
* openSUSE 11.3
* Ubuntu (10.04/Lucid and 10.10/Maverick)

A number of other platforms are supported through user contributions.
It is easy to extend support to other similar platforms. Please let us
know if you are interested, and we can work together towards supporting
your platform.

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
`FFC <http://launchpad.net/ffc>`_,
`FIAT <http://launchpad.net/fiat>`_,
`Instant <http://launchpad.net/instant>`_,
`Viper <http://launchpad.net/fenics-viper>`_, and
`UFL <http://launchpad.net/ufl>`_.
You will also need to install the C++/Python packages
`UFC <http://launchpad.net/ufc>`_ and
`DOLFIN <http://launchpad.net/dolfin>`_.
Instructions for how to install these packages are given below.

Python packages (FFC, FIAT, Instant, Viper and UFL)
===================================================

#. First, download the source code for FFC, FIAT, Instant, Viper and
   UFL. Then for each of these packages, follow the instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following command::

    sudo python setup.py install

This will install the packages system wide on your computer. You may
easily change the installation directory. For example, if you do not
have super-user access, you may wish to install the packages in a
directory named ``local`` in your home directory. This may be
accomplished by running the following command::

    python setup.py install --prefix=~/local

C++/Python packages (DOLFIN and UFC)
====================================================

Building and installing DOLFIN and UFC require CMake version >= 2.8.

#. First, download the source code for UFC and DOLFIN. Then for each of
   these packages, follow the instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following commands::

    cmake .
    make
    sudo make install

This will install the packages system wide on your computer. You may easily
change the installation directory. For example, you may wish to install
the packages in a directory named ``local`` in your home directory. This
may be accomplished by replacing the first of the above commands by::

    cmake -DCMAKE_INSTALL_PREFIX=~/local .

It is also possible (and usually recommended) to build DOLFIN 'out of
source'.  In the directory where you wish to build DOLFIN, the build can
be configured by::

    cmake -DCMAKE_INSTALL_PREFIX=<prefix> <source_path>

where <source_path> is the path to the DOLFIN source.
You may also want to consider using a graphical front end for CMake such
as either ``cmake-gui`` or ``ccmake``. These both provide a simple way to
configure the installation.

During the configuration phase of DOLFIN (calling ``cmake``, ``cmake-gui``, or
``ccmake``), you will be notified of any missing dependencies. If a required
package is missing, you will need to install that package and configure DOLFIN
again. If an optional package is missing, you may choose to continue with the
installation but some functionality may be missing. The build system will list
both found and missing optional dependencies at the end of the configuration
process.

***************
Troubleshooting
***************

If for some reason you were unable to obtain, build, or install a FEniCS
component, please let us know. The best place to start is the `Launchpad
Answers <https://help.launchpad.net/Answers>`_ page for the project that
causes the problem:

* `DOLFIN <http://answers.launchpad.net/dolfin>`__
* `Dorsal <https://answers.launchpad.net/dorsal>`__
* `FFC <https://answers.launchpad.net/ffc>`__
* `FIAT <https://answers.launchpad.net/fiat>`__
* `Instant <https://answers.launchpad.net/instant>`__
* `Viper <https://answers.launchpad.net/fenics-viper>`__
* `UFC <https://answers.launchpad.net/ufc>`__
* `UFL <https://answers.launchpad.net/ufl>`__

If all else fails, send a message to the fenics@lists.launchpad.net
mailing list.

If you have suggestions for improving the information found on this page
with respect to the installation procedure of FEniCS, you are more than
welcome to file a bug report or register a blueprint on `FEniCS Documentation
<https://launchpad.net/fenics-doc>`_.
