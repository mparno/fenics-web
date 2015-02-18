.. _installation_using_conda:

###############################
Binary installation using conda
###############################

_Anaconda is a Python distribution focused on scientific computing from the
company Continuum Analytics. It is centered around conda_, an open source,
language-agnostic package manager. Its distinct features are:

* It downloads binary packages in all platforms.
* It is built around environments, which are first-class citizens and way more powerful than virtualenvs.
* It is language-agnostic - you can use it to install things like boost or vtk too.

FEniCS packages for conda are user-contributed and support any Linux 64-bit
distribution. As they are not part of the official
Anaconda packages, so you have to use a binstar_ channel to install them:

    conda create --name fenics27 python=2.7
    source activate fenics27
    conda install fenics --channel juanlu001

This will download all FEniCS dependencies and install them quickly.
Alternatively, you can add this channel to your conda settings so it will search
for packages there too:

    conda config --add channels juanlu001
    source activate fenics27
    conda install fenics  # This will work now

If you happen to be using some old system, FEniCS 1.4.0 packages are provided
compiled in CentOS 6. These can be installed from another binstar channel:

    conda install "fenics=1.4.0" mkl --channel https://conda.binstar.org/juanlu001/channel/fenics:1.4.0:centos

Or alternatively:

    conda config --add channels https://conda.binstar.org/juanlu001/channel/fenics:1.4.0:centos
    conda install "fenics=1.4.0"

.. note::

    FEniCS binaries from this channel are linked to Continuum's Intel MKL, which
    is available as a paid-for product or under an `academic license`_. If you want
    to use your own MKL or BLAS/LAPACK or simply don't want to use the default MKL
    you will have to build your own packages, see below.

.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _conda: http://conda.io/
.. _binstar: https://binstar.org/

.. _`academic license`: https://store.continuum.io/cshop/academicanaconda

############################################
Building custom conda packages using recipes
############################################

conda is not only a binary package manager, but also includes a build system
called `conda-build`_ which in turn makes it straightforward to build your
own conda packages, either for private use or for the general public using
binstar channels. A conda package is created using a *conda recipe* as explained
in the documentation. Such recipes for FEniCS can be found in the
`fenics-recipes repository <https://github.com/juanlu001/fenics-recipes>`_.
As explained there, to recreate this packages in your system from source
you can use `conda-build`:

    conda install conda-build
    conda build boost eigen3 petsc petsc4py instant ufl fiat ffc dolfin fenics --python 27
    conda install fenics mkl --use-local

.. _`conda-build`: http://conda.pydata.org/docs/build.html

###############
Troubleshooting
###############

"Error: No packages found matching ..." when installing
-------------------------------------------------------

You probably forgot to specify the channel, or to add them to your conda
configuration, as explained above.

"ImportError: /lib64/libc.so.6: version `GLIBC_2.14' not found"
---------------------------------------------------------------

Maybe you installed the latest FEniCS version but your system is too old.
Try to install FEniCS 1.4.0 from the channel
https://conda.binstar.org/juanlu001/channel/fenics:1.4.0:centos as explained
above.

"ImportError: ... cannot open shared object file: No such file or directory"
----------------------------------------------------------------------------

There is some sort of linking problem in your system. Perhaps you have
to install some of the Boost system requirements, specially libbz2. To diagnose
this problem you can use the `ldd` utility:

    ldd <INSTALL_PATH>/anaconda/envs/fenics27/lib/python2.7/site-packages/dolfin/cpp/_common.so

My error is not listed above
----------------------------

Please open an issue at https://github.com/juanlu001/fenics-recipes/issues.
In the meanwhile, you can try and build your own conda packages from the recipes.
