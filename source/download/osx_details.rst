.. _osx_details:

########
Download
########

This release includes DOLFIN 1.3.0, FFC 1.3.0, FIAT 1.3.0, Instant 1.3.0,
UFC 2.3.0 and UFL 1.3.0.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |mac-icon| FEniCS 1.3.0 (Mac OS X 10.7 binary)

      - `fenics-1.3.0-osx10.7.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.7.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.7.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.7.dmg.asc>`__)

    * - |mac-icon| FEniCS 1.3.0 (Mac OS X 10.8 binary)

      - `fenics-1.3.0-osx10.8.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.8.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.8.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.8.dmg.asc>`__)

    * - |mac-icon| FEniCS 1.3.0 (Mac OS X 10.9 binary)

      - `fenics-1.3.0-osx10.9.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.9.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.9.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.3.0-osx10.9.dmg.asc>`__)

############
Requirements
############

* 64-bit capable Intel processor.
* Mac OS X 10.7 Lion, 10.8 Mountain Lion or 10.9 Mavericks.
* An up-to-date version of the Xcode Command Line Developer Tools. On OS
  X 10.7 and 10.8, install Xcode from App Store and select the Xcode
  Command Line Tools from the Downloads tab within the Xcode Preferences
  menu. On OS X 10.9, simply run ``xcode-select --install`` and click on
  the Install button.
* The `X11 Libraries <http://xquartz.macosforge.org/trac/wiki>`__ must
  be installed on OS X 10.8.
* The FEniCS binary is built against system Python and will *not work*
  with Python from `python.org <http://python.org>`__, Python from
  MacPorts or similar.

############
Installation
############

Click on the link above to download the binary package for your version
of OS X. The installation should be as simple as dragging the FEniCS
icon into the Applications folder. When the installation is complete,
there are two ways to use this binary:

* Click on the FEniCS icon in the Applications folder. This will bring
  up a terminal with everything set up to work with FEniCS.

* Add the following line to the ``.profile`` file in your home
  directory::

    source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf

  This will make FEniCS available whenever you start a new terminal.

##############
Older releases
##############

Older versions of the Mac OS X binary are also :ref:`available
<older_releases>`.
