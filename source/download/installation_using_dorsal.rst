.. _installation_using_dorsal:

###################################
Automated installation using Dorsal
###################################

If you are not running one of the operating systems we provide binary
packages for, you need to download and compile FEniCS and all its
dependencies. Luckily, this can be handled easily using `Dorsal
<https://bitbucket.org/fenics-project/dorsal>`__.  Dorsal is a simple
shell script that automates the process of fetching, compiling, and
installing the various FEniCS sub-projects and their requisite
dependencies on Unix-like operating systems.

The following platforms are officially supported by Dorsal:

* Debian GNU/Linux 6.0, Unstable
* Fedora 13, 14, 15
* Gentoo Linux
* Mac OS X 10.6, 10.7 (with MacPorts)
* openSUSE 11.3, 11.4
* Red Hat Enterprise Linux 6
* Ubuntu 12.04 LTS, 12.10, 13.04

A number of other platforms are supported through user contributions.
It is easy to extend support to other similar platforms. Please let us
know if you are interested, and we can work together towards supporting
your platform.

In order to install FEniCS using Dorsal, simply follow these steps:

#. Fetch the most recent snapshot of Dorsal from
   `Bitbucket <https://bitbucket.org/fenics-project/dorsal>`_.
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
