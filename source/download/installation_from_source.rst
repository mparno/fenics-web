.. _installation_from_source:

*******************************
Manual installation from source
*******************************

You can also build and install FEniCS components manually from the source code.
You will need to install the Python packages
`FFC <https://bitbucket.org/fenics-project/ffc>`_,
`FIAT <https://bitbucket.org/fenics-project/fiat>`_,
`Instant <https://bitbucket.org/fenics-project/instant>`_,
`UFL <https://bitbucket.org/fenics-project/ufl>`_.
You will also need to install the C++/Python packages
`UFC <https://bitbucket.org/fenics-project/ufc-deprecated>`_ and
`DOLFIN <https://bitbucket.org/fenics-project/dolfin>`_.
Instructions for how to install these packages are given below.

Python packages (FFC, FIAT, Instant and UFL)
============================================

#. First, download the source code for `FFC
   <https://bitbucket.org/fenics-project/ffc>`_, `FIAT
   <https://bitbucket.org/fenics-project/fiat>`_, `Instant
   <https://bitbucket.org/fenics-project/instant>`_ and `UFL
   <https://bitbucket.org/fenics-project/ufl>`_. Then for each of these
   packages, follow the instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following command::

    sudo python setup.py install

This will install the packages system wide on your computer. You may
easily change the installation directory. For example, if you do not
have super-user access, you may wish to install the packages in a
directory named ``local`` in your home directory. This may be
accomplished by running the following command::

    python setup.py install --prefix=~/local

C++/Python packages (DOLFIN and UFC)
====================================

Building and installing DOLFIN and UFC require CMake version >= 2.8.

#. First, download the source code for `UFC
   <https://bitbucket.org/fenics-project/ufc-deprecated>`_ and `DOLFIN
   <https://bitbucket.org/fenics-project/dolfin>`_. Then for each of
   these packages, follow the instructions below.
#. Uncompress the archive to a convenient location.
#. Navigate to this folder and run the following commands::

    cmake .
    make
    sudo make install

This will install the packages system wide on your computer. You may easily
change the installation directory. For example, you may wish to install
the packages in a directory named ``local`` in your home directory. This
may be accomplished by replacing the first of the above commands by::

    cmake -DCMAKE_INSTALL_PREFIX=~/local .

It is also possible (and usually recommended) to build DOLFIN 'out of
source'.  In the directory where you wish to build DOLFIN, the build can
be configured by::

    cmake -DCMAKE_INSTALL_PREFIX=<prefix> <source_path>

where <source_path> is the path to the DOLFIN source.
You may also want to consider using a graphical front end for CMake such
as either ``cmake-gui`` or ``ccmake``. These both provide a simple way to
configure the installation.

During the configuration phase of DOLFIN (calling ``cmake``, ``cmake-gui``, or
``ccmake``), you will be notified of any missing dependencies. If a required
package is missing, you will need to install that package and configure DOLFIN
again. If an optional package is missing, you may choose to continue with the
installation but some functionality may be missing. The build system will list
both found and missing optional dependencies at the end of the configuration
process.
