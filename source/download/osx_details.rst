.. _osx_details:

######################################
Installation instructions for Mac OS X
######################################

This release includes DOLFIN 1.5.0, FFC 1.5.0, FIAT 1.5.0, Instant 1.5.0
and UFL 1.5.0.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |mac-icon| FEniCS 1.5.0 (Mac OS X 10.10 binary)

      - `fenics-1.5.0-p2-osx10.10.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.10.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.10.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.10.dmg.asc>`__)

    * - |mac-icon| FEniCS 1.5.0 (Mac OS X 10.9 binary)

      - `fenics-1.5.0-p2-osx10.9.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.9.dmg>`__
        (`md5
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.9.dmg.md5>`__,
        `sig
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.5.0-p2-osx10.9.dmg.asc>`__)

************
Requirements
************

* Mac OS X 10.9 Mavericks or 10.10 Yosemite.
* An up-to-date version of the Xcode Command Line Developer Tools. To
  install, simply run ``xcode-select --install`` in a terminal and click
  on the Install button.
* The FEniCS binary is built against system Python and will *not work*
  with Python from `python.org <http://python.org>`__, Python from
  MacPorts or similar.

************
Installation
************

The installation should be as simple as dragging the FEniCS icon into
the Applications folder. When the installation is complete, there are
two ways to use this binary:

* Click on the FEniCS icon in the Applications folder. This will bring
  up a terminal with everything set up to work with FEniCS.

* Add the following line to the ``.profile`` file in your home
  directory::

    source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf

  This will make FEniCS available whenever you start a new terminal.

**************
Older releases
**************

Older versions of the Mac OS X binary are also :ref:`available
<older_releases>`.
