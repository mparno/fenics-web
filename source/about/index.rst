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
are documented in a number of :doc:`research articles
<../citing/index>` and a :doc:`book <../book/index>`. You can also
learn more about FEniCS from the following `presentations
<http://fenicsproject.org/pub/presentations/>`__.

********
Features
********

FEniCS has an :ref:`extensive list of features <features>` for
automated, efficient solution of differential equations, including
automated solution of variational problems, automated error control
and adaptivity, a comprehensive library of finite elements, high
performance linear algebra and many more.

**********
Components
**********

FEniCS is organized as a collection of interoperable :ref:`components
<about_components>` that together form the FEniCS Project. These
components include the problem-solving environment :ref:`DOLFIN
<about_components_dolfin>`, the form compiler :ref:`FFC
<about_components_ffc>`, the finite element tabulator :ref:`FIAT
<about_components_fiat>`, the just-in-time compiler :ref:`Instant
<about_components_instant>`, the code generation interface :ref:`UFC
<about_components_ufc>`, the form language :ref:`UFL
<about_components_ufl>` and a range of :ref:`additional components
<about_components_additional>`.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <img src="../_images/fenics-map.png" class="img-responsive" alt="Fenics Map">
    </div>
  </div>


Building on these components, software specialized to solving
different problems are organised into separate :doc:`applications
<../applications/index>`.

************
Contributors
************

The FEniCS Project is developed by :ref:`researchers from a number of
research institutes <team>` from around the world. The following
research institutes (in alphabetical order) contribute to the FEniCS
Project:

* `Baylor University <http://www.baylor.edu>`__
* `Chalmers University of Technology <http://www.chalmers.se>`__
* `KTH Royal Institute of Technology <http://www.kth.se>`__
* `Simula Research Laboratory <http://www.simula.no>`__
* `University of Cambridge <http://www.cam.ac.uk>`__
* `University of Chicago <http://www.uchicago.edu/index.shtml>`__

Contributions have also been made by researchers from `Delft
University of Technology <http://www.tudelft.nl>`__, `Argonne National
Laboratory <http://www.anl.gov>`__ and many other research institutes.
A full list of contributors is maintained as part of the source code
of each FEniCS component.

The following video illustrates the development of the FEniCS Project
since its inception in 2003.

.. raw:: html

  <div class="container-fluid">
    <div class="row">
      <div class="center-block" style="max-width:500px">
        <div class="embed-responsive embed-responsive-16by9" style=>
          <iframe class="embed-responsive-item" width="50%" src="http://www.youtube.com/embed/0E6DGUbRrO4"></iframe>
        </div>
      </div>
    </div>
  </div>

.. _about_license:

*******
License
*******

All FEniCS :ref:`core components <about_components_core>` are licensed
under the `GNU LGPL <http://www.gnu.org/licenses/lgpl.html>`__ as
published by the `Free Software Foundation <http://www.fsf.org>`__,
either version 3 of the license, or (at your option) any later
version.

All other FEniCS components are licensed under either the `GNU GPL
<http://www.gnu.org/licenses/gpl.html>`__ or the `GNU LGPL
<http://www.gnu.org/licenses/lgpl.html>`__, either version 3 of the
license, or (at your option) any later version.

Authors and institutions have given their consent to licensing under
`these terms <http://www.fenicsproject.org/pub/copyright>`__.

*****************
About these pages
*****************

These pages have been created by the `FEniCS Team
<https://bitbucket.org/fenics-project/profile/members>`__ with the
help of `Mattias Schläger <http://www.sch-form.com>`__ who designed
the `graphical profile
<http://www.fenicsproject.org/pub/graphics>`__. The pages are built
using the `Sphinx documentation system <http://sphinx.pocoo.org>`__ in
combination with some homebrew scripting (for extraction of C++
documentation in particular). The sources for these pages are
maintained on `Bitbucket
<https://bitbucket.org/fenics-project/fenics-web>`__. Comments and bug
reports are welcome as always. If you find something is in error or
missing, `please file a bug report
<https://bitbucket.org/fenics-project/fenics-web/issues>`__ on
Bitbucket.

.. toctree::
   :hidden:

   components
   features
   team
