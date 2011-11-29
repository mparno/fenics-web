.. _osx_details:

########
Download
########

This release includes DOLFIN 1.0-rc2, FFC 1.0-rc1, FIAT 1.0-beta,
Instant 1.0-beta, UFC 2.0.4, UFL 1.0-rc1, and Viper 1.0-beta.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |mac-icon| FEniCS 1.0-rc2 (Mac OS X 10.6 binary)

      - `fenics-1.0-rc2-osx10.6.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.0-rc2-osx10.6.dmg>`__

    * - |mac-icon| FEniCS 1.0-rc2 (Mac OS X 10.7 binary)

      - `fenics-1.0-rc2-osx10.7.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-1.0-rc2-osx10.7.dmg>`__

    * - |archive-icon| FEniCS 1.0-rc2 (Sources for Mac OS X 10.6)

      - `fenics-source-1.0-rc2-osx10.6.tar.gz
        <http://www.fenicsproject.org/pub/software/fenics/fenics-source-1.0-rc2-osx10.6.tar.gz>`__

    * - |archive-icon| FEniCS 1.0-rc2 (Sources for Mac OS X 10.7)

      - `fenics-source-1.0-rc2-osx10.7.tar.gz
        <http://www.fenicsproject.org/pub/software/fenics/fenics-source-1.0-rc2-osx10.7.tar.gz>`__

#########################
Installation instructions
#########################

The FEniCS binary for Mac OS X runs on Intel 10.6 and 10.7 only. This
means that OS X 10.5 and older and the PowerPC architecture *are not
supported*. Moreover, this package is built against system Python and
will *not work* with MacPorts Python or similar. XCode (available from
App Store) is required to run the binary.

Click on the link above to download the binary package. The installation
should be as simple as dragging the FEniCS icon into the Applications
folder. When the installation is complete, there are two ways to use
this binary:

* Click on the FEniCS icon in the Applications folder. This will bring
  up a terminal with everything set up to work with FEniCS.

* Add the following line to the ``.profile`` file in your home
  directory::

    source /Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf

  This will make FEniCS available whenever you start a new terminal.

**************
Older releases
**************

Older versions of the Mac OS X binary are available :ref:`here
<older_releases>`.
