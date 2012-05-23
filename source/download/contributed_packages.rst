.. _contributed_packages:

####################
Contributed packages
####################

The packages listed on this page are user contributed packages and not
provided by the FEniCS team. If you have problems with them, please
contact the package maintainer first.

If you have packages or information about other systems, it is always
appreciated if you let us know. Please send the information to the
`fenics-web mailing list <fenics-web@lists.launchpad.net>`__, or even
better, send us a patch or make a merge request against `lp:fenics-web
<https://code.launchpad.net/~fenics-web-core/fenics-web/main>`__.

**********
Arch Linux
**********

FEniCS is available in the `Arch User Repository
<https://wiki.archlinux.org/index.php/Arch_User_Repository>`__ (AUR),
thanks to Myles English. The packages follows the trunk versions of the
FEniCS projects.

Start by running the following command to install the required packages
from the official repository::

    sudo pacman -S bzr cmake swig boost libxml2 lapack blas python2-numpy python2-scientificpython

Then build and install all these packages from the AUR in the following
order (along with any other dependencies they may need): ``parmetis``,
``metis``, ``petsc``, ``slepc``, ``trilinos``, ``instant-bzr``,
``fiat-bzr``, ``ufl-bzr``, ``ufc-bzr``, ``ffc-bzr``, ``viper-bzr``,
``dolfin-bzr``.
