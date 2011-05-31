.. Developer information.

.. _developers:

#####################
Developer information
#####################

This page contains information for FEniCS developers, including an
overview of the organization of the FEniCS project, how to use
Launchpad and Bazaar, and how to write code and documentation.

************
Organization
************

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. Each component is developed by one
or more authors. This means that each component can be developed at
its own pace, but we strive to make periodic and coordinated releases
of all components to ensure interoperability between the components.

Initially, FEniCS consisted of just two components (DOLFIN and FIAT)
but over time, several new components have been added and FEniCS now
consists of more than 10 individual components. Some of these
components (such as FIAT and UFC) have matured and reached a more
stable state, while others are changing at a faster pace. Currently,
most development takes place in DOLFIN, the C++ and Python interface
of FEniCS.

.. toctree::
   :numbered:
   :maxdepth: 2

   using_launchpad
   using_bazaar
   contributing_code
   automated_testing
   writing_documentation
