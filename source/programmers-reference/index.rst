.. The programmer's reference for FEniCS goes here.

.. _programmers_reference_index:

#############################
FEniCS Programmer's Reference
#############################

The FEniCS Programmer's Reference documents the details of the FEniCS
user interface. The FEniCS user interface is available in two
different flavors: both as a Python module and as a C++ library.

Both interfaces are implemented as part of `DOLFIN
<http://www.fenics.org/dolfin>`_.

Python interface
================

To use FEniCS from Python, users need to import functionality from the
DOLFIN Python module. In the simplest caset, one includes all
functionality from the Python module named ``dolfin``:

.. code-block:: python

    from dolfin import *

For documentation of the DOLFIN Python module, see the `Python Programmer's Reference <python/index.html>`_.

C++ interface
=============

To use FEniCS from C++, users need to include one or more header files
from the DOLFIN C++ library. In the simplest case, one includes the
header file ``dolfin.h``, which in turn includes all other DOLFIN
header files:

.. code-block:: c++

    #include <dolfin.h>

    using namespace dolfin;

    int main()
    {

      return 0;
    }

For documentation of the DOLFIN C++ library, see the `C++ Programmer's Reference <cpp/index.html>`_.
