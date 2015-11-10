.. _contributed_packages:

####################
Contributed packages
####################

The packages listed on this page are user contributed packages and not
provided by the FEniCS team. If you have problems with them, please
contact the package maintainer first.

If you have packages or information about other systems, it is always
appreciated if you let us know. Please send the information to the
`fenics mailing list <fenics-dev@googlegroups.com>`__, or even better,
send us a patch or make a `pull request
<https://bitbucket.org/fenics-project/fenics-web>`__.

**********
Arch Linux
**********

FEniCS is available in the `Arch User Repository
<https://wiki.archlinux.org/index.php/Arch_User_Repository>`__ (AUR),
thanks to Myles English. The packages follow the trunk versions of the
FEniCS projects.

Start by running the following command to install the required packages
from the official repository::

    sudo pacman -S bzr cmake swig boost libxml2 lapack blas python2-numpy

Then build and install all these packages from the AUR in the following
order (along with any other dependencies they may need):
``python2-scientificpython``, ``metis``, ``parmetis``, ``petsc``,
``slepc``, ``trilinos``, ``instant-git``, ``fiat-git``, ``ufl-git``,
``ffc-git``, ``dolfin-git``.

********
openSUSE
********

FEniCS is available in the `openSUSE Science Repository
<http://download.opensuse.org/repositories/science/>`__, thanks to
Sebastien Corot.

First we need to add the Science Repository if this has not been already
added. Start by opening YaST2 and select *Software Repositories*. Then
click on the button labeled *Add*. This brings up a screen labeled
*Media Type*. From there select ``HTTP`` and click on *Next*. This brings
up a screen labeled *Repository URL*. For the Science Repository, fill
in::

  Repository Name: openSUSE Science
  URL of the Repository: http://download.opensuse.org/repositories/science/openSUSE_12.3/

Make sure to replace the version number with the correct version of
openSUSE. Then click on *Next*. That should bring up a screen that is
labeled *Adding New Repository*, and then return you to the list of
configured repositories for your system. This list should now contain
the openSUSE Science Repository.

Now, to install DOLFIN, one goes into YaST2 and either searches for
``dolfin`` or selects ``Repositories`` and then selects ``dolfin`` from
the list of available packages. YaST2 should also take care of selecting
any dependencies.
