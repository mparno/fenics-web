.. General introduction to the FEniCS Project.

.. _About:

########################
About the FEniCS Project
########################

The FEniCS Project is a collaborative project for the development of
innovative concepts and tools for automated scientific computing, with
a particular focus on automated solution of differential equations by
finite element methods.

The methodology and software developed as part of the FEniCS Project
are documented in a number of :doc:`research articles <../citing/index>`
and a :doc:`book <../book/index>`.

************
Contributors
************

The FEniCS Project is developed by researchers from a number of
research institutes from around the world. The following research
institutes contribute significant resources to the FEniCS Project:

* `Simula Research Laboratory <http://www.simula.no>`__
* `University of Cambridge <http://www.cam.ac.uk>`__
* `University of Chicago <http://www.uchicago.edu/index.shtml>`__
* `Texas Tech University <http://www.ttu.edu/>`__
* `KTH Royal Institute of Technology <http://www.kth.se>`__

Contributions have also been made by researchers from `Chalmers
University of Technology <http://www.chalmers.se>`__, `Delft
University of Technology <http://www.tudelft.nl>`__, `Argonne National
Laboratory <http://www.anl.gov>`__ and many other research institutes.
A full list of contributors is maintained as part of the source code
of each FEniCS component.

.. note::
    Johannes' revision-control-history videos can be featured here.

**********
Components
**********

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. These components include the
problem-solving environment :ref:`DOLFIN <about_components_dolfin>`,
the form compiler :ref:`FFC <about_components_ffc>`, the finite
element tabulator :ref:`FIAT <about_components_fiat>`, the
just-in-time compiler :ref:`Instant <about_components_instant>`, the
code generation interface :ref:`UFC <about_components_ufc>`, the form
lanuage :ref:`UFL <about_components_ufl>` and a range of
:ref:`additional components <about_components_additional>` and
:ref:`applications <about_components_applications>`.

.. note::
    I think applications should be removed from here and featured as
    part of their :ref:`own page <apps>`.

An full overview of the list of FEniCS components is presented
:ref:`here <about_components>`.

.. _about_license:

*******
License
*******

All FEniCS :ref:`core components <about_components_core>` are licensed
on the `GNU LGPL <http://www.gnu.org/licenses/lgpl.html>`__ as
published by the `Free Software Foundation <http://www.fsf.org>`__,
either version 3 of the license, or (at your option) any later
version.

All other FEniCS components are licensed under either the `GNU GPL
<http://www.gnu.org/licenses/gpl.html>`__ or the `GNU LGPL
<http://www.gnu.org/licenses/lgpl.html>`__, either version 3 of the
license, or (at your option) any later version.

Authors and institutions have given their consent to licensing under
these terms `here <http://www.fenicsproject.org/pub/copyright>`__.

.. toctree::
   :hidden:

   components

.. note::
    In which order should we list the research institutes? In order of
    contribution, or chronological?
