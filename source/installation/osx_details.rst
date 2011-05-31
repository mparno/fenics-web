########
Download
########

This release includes DOLFIN 0.9.11, FFC 0.9.10, FIAT 0.9.9, Instant
0.9.10, UFC 2.0.1, UFL 0.9.1, and Viper 0.4.7.

.. include:: icons.rst

.. tabularcolumns:: |c|c|

.. list-table::
    :widths: 50, 50
    :header-rows: 0
    :class: center

    * - |mac-icon| FEniCS 11.05 (Mac OS X binary)

      - `fenics-11.05-osx10.6.dmg
        <http://www.fenicsproject.org/pub/software/fenics/fenics-11.05-osx10.6.dmg>`__

    * - |archive-icon| FEniCS 11.05 (Sources for Mac OS X)

      - `fenics-source-11.05-osx.tar.gz
        <http://www.fenicsproject.org/pub/software/fenics/fenics-source-11.05-osx.tar.gz>`__

#########################
Installation instructions
#########################

The FEniCS binary for Mac OS X runs on Intel 10.6 only. This means that
10.5 and older and the PowerPC architecture *are not
supported*. Moreover, this package is built against system Python and
will *not work* with MacPorts Python or similar. XCode 3 (available on
the OS X install disc) is required to run the binary.

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

Older versions of the Mac OS X binary are available `here
<older_releases.html>`__.
