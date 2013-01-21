.. _osx_details:

########
Download
########

This release includes DOLFIN 1.1.0, FFC 1.1.0, FIAT 1.1, Instant 1.1.0,
UFC 2.1.0 and UFL 1.1.0.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |mac-icon| FEniCS 1.1.0 (Mac OS X 10.6 binary)

      - `fenics-1.1.0-osx10.6.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.6.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.6.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.6.dmg.asc>`__)

    * - |mac-icon| FEniCS 1.1.0 (Mac OS X 10.7 / 10.8 binary)

      - `fenics-1.1.0-osx10.7.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.7.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.7.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.1.0-osx10.7.dmg.asc>`__)

############
Requirements
############

* 64-bit capable Intel processor.
* Mac OS X 10.6 or greater.
* An up-to-date version of Xcode. This can be installed from the OS X
  install disc for Snow Leopard users, or from App Store for Lion or
  Mountain Lion users.
* The Xcode Command Line Tools must be installed from the Downloads tab
  within the Xcode Preferences menu for Lion or Mountain Lion users.
* The `X11 Libraries <http://xquartz.macosforge.org/trac/wiki>`__ must
  be installed on Mountain Lion.
* The FEniCS binary is built against system Python and will *not work*
  with Python from `python.org <http://python.org>`__, MacPorts Python
  or similar.

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
