.. _windows_details:

#####################################
Installation instructions for Windows
#####################################

.. warning::

    The FEniCS Windows binary is not up-to-date. To install the latest
    FEniCS release on Windows, we currently recommend to install `Ubuntu
    <http://www.ubuntu.com>`__ in a virtual machine (e.g. using
    `VirtualBox <http://www.virtualbox.org>`__) and then install FEniCS
    from the :ref:`Ubuntu PPA <ubuntu_ppa>`.

The latest (and outdated) release for Windows includes DOLFIN 1.0.0, FFC 1.0.0, FIAT 1.0.0,
Instant 1.0.0, UFC 2.0.5, UFL 1.0.0, and Viper 1.0.0.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |windows-icon| FEniCS 1.0.0 (Windows installer)

      - `fenics-1.0.0-mingw32.exe
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.0.0-mingw32.exe>`__

    * - |archive-icon| FEniCS 1.0.0 (Sources for Windows)

      - `fenics-source-1.0.0-win.zip
        <http://www.fenicsproject.org/pub/software/fenics/fenics-source-1.0.0-win.zip>`__

The Windows installer will install everything needed to run FEniCS on
Windows, including MinGW compilers, Python, CMake, SWIG, and others. It
does not currently include PETSc or Trilinos, which makes it unsuitable
for solving any real problems. The installer has been verified to run on
XP, Vista, and on Windows 7.

Click the link above to download the Windows installer. The installation
procedure is simple, just double-click on the file and follow the
instructions. When the installation is complete, simply start the FEniCS
command shell from the start menu and you are ready to start working.

**************
Older releases
**************

Older versions of the Windows installer are also :ref:`available
<older_releases>`.
