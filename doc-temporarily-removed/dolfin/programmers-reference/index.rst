.. The programmer's reference for DOLFIN goes here.

.. _programmers_reference_index:

#############################
DOLFIN programmer's reference
#############################

The FEniCS programmer's reference documents the FEniCS
interface. The FEniCS interface is available in two
different flavors: as a C++ library and as a Python module.
Both interfaces are implemented as part of `DOLFIN
<http://www.fenics.org/dolfin>`_.

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

For documentation of the DOLFIN C++ library, see the
:ref:`programmers_reference_cpp_index`.

Python interface
================

To use FEniCS from Python, users need to import functionality from the
DOLFIN Python module. In the simplest case, one includes all
functionality from the Python module named ``dolfin``:

.. code-block:: python

    from dolfin import *

For documentation of the DOLFIN Python module, see the
:ref:`doc_dolfin_programmers_reference_python_index`.

.. toctree::
    :hidden:

    cpp/index
    python/index

