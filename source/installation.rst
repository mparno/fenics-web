.. How to download and install FEniCS projects.

.. _installation:

############
Installation
############

Some general remarks on installing FEniCS software. Where to find it, how to
get it and links to Launchpad.

************************************
Installation for Debian/Ubuntu users
************************************

Explain how to install the packages.

*************************
Installation using Dorsal
*************************

The easiest way to install FEniCS on a UNIX-like operating environment
is to use Dorsal. Dorsal is a simple shell script that automates the
process of fetching, compiling and installing the various FEniCS
sub-projects and their requisite dependencies. It currently supports
the following platforms:

* Debian GNU/Linux 4.0, 5.0, unstable
* Fedora 10, 11, 12          
* Mac OS X 10.4, 10.5, 10.6
* openSUSE 11.1, 11.2      
* Red Hat Enterprise Linux 4, 5 
* Ubuntu 8.10, 9.04, 9.10, 10.04 LTS
* Gentoo (and Sabayon) Linux

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
and their corresponding dependencies.

Don't forget to follow any post-build instructions before rushing off to try
the demos!


***
FFC
***

Install FFC

******
DOLFIN
******

Install DOLFIN

***
UFL
***

Install UFL

***
UFC
***

Install UFC

