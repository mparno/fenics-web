.. Developer information.

.. _developers:

############
Contributing
############

The FEniCS Project is an open source project where everyone is invited
to contribute.

The project is organized as :ref:`a collection of interoperable
components <About>` that together form the FEniCS Project. Each
component is developed by one or more authors. This means that each
component can be developed at its own pace. At the same time, we
strive to make periodic and coordinated releases of all components to
ensure interoperability between the components.

* :ref:`How to obtain the FEniCS source code <obtaining_the_source_code>`
* :ref:`How to contribute code to FEniCS <contributing_code>`
* :ref:`How to contribute a featured FEniCS story <contributing_featured_fenics_stories>`

.. _obtaining_the_source_code:

*************************
Obtaining the source code
*************************

The FEniCS Project uses `Bitbucket <https://bitbucket.org>`_ as its
main development platform. Each FEniCS component has a `Git
<http://git-scm.com/>`_ repository on Bitbucket that contains all
source code (including the entire development history). All FEniCS
projects are collected under the `FEniCS Project on Bitbucket
<https://bitbucket.org/fenics-project>`_.

The repositories are readable for everyone, but write access is only
granted to the members of the core teams.

======================================
Accessing the development repositories
======================================

To access the development repositories, you need the revision control
system `Git <http://git-scm.com/>`__.

You can access the development repository of any FEniCS project via
the ``git`` terminal command. For example, to check out the source
code for DOLFIN, type::

    git clone https://bitbucket.org/fenics-project/dolfin.git

or if you have an ssh key associated with your Bitbucket account::

    git clone git@bitbucket.org:fenics-project/dolfin.git

========================
Notifications of updates
========================

Developers should subscribe to notifications of changes made to the
source code by visiting each of `the FEniCS repositories on Bitbucket
<https://bitbucket.org/fenics-project>`__ and clicking the 'follow'
button.

======================
FEniCS Developer Tools
======================

Developers should take a look at the `FEniCS Developer Tools
<https://bitbucket.org/fenics-project/fenics-developer-tools>`__
repository. This contains scripts that are very useful for developers
in setting up a good development environment for FEniCS. In
particular, consider using the scripts ``fenics-install-all.sh`` and
``fenics-install-component.sh``.

.. _contributing_code:

*****************
Contributing code
*****************

Unless you are a core developer, you will not be able to upload any
code to the main code repositories directly. Instead you will have to
create your own branch and make a merge request on Bitbucket. If the
code is accepted, the patch or branch will be merged into the main
branch by a core developer.

If you are not familiar with `Git <http://git-scm.com/>`__, the
distributed revision control system used for all FEniCS components, a
good starting point is `Bitbucket 101
<https://confluence.atlassian.com/display/BITBUCKET/Bitbucket+101>`__
or the `Git user documentation <http://git-scm.com/documentationl>`__.

============
Git workflow
============

FEniCS development follows the `gitworkflows
<https://www.kernel.org/pub/software/scm/git/docs/gitworkflows.html>`__
model (with the exception of 'pu' branches). Developers should read up
on the gitworkflows model and understand the role of 'master', 'next'
and 'topic branches'. The same workflow is used by the developers of
`PETSc <http://www.mcs.anl.gov/petsc/>`__. The `PETSc Wiki
<https://bitbucket.org/petsc/petsc/wiki/Home>`__ has some good
information on both Git usage and the gitworkflows model.

A summary of useful Git commands for some common use cases can be
found in the `Git cookbook for FEniCS developers
<https://bitbucket.org/fenics-project/dolfin/wiki/Git%20cookbook%20for%20FEniCS%20developers>`__.

=======
Testing
=======

Before submitting any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying. If you are
submitting code for a new feature, it is expected that the new feature
is accompanied by a suitable set of unit tests.

==================
FEniCS style guide
==================

We ask all developers and contributors to adhere to a common style
guide. This makes the job easier for maintainers who need to review,
edit and maintain the FEniCS code base. The FEniCS coding style for
Python code adheres to the `PEP-8 style guide
<http://www.python.org/dev/peps/pep-0008/>`_ although it is not
strictly enforced. Guidelines for writing C++ code are given in the
DOLFIN doc/ directory.

..
   Patches that don't follow the correct coding style will likely be
   rejected, as the maintainer responsible for reviewing the patch must
   otherwise make additional efforts to edit the patch to follow the
   coding style.

===============================
Copyright and licensing consent
===============================

Before your contribution can be accepted into FEniCS, you must sign a
`copyright consent form
<http://fenicsproject.org/pub/copyright/forms/>`_.  Ideally, both you
and your employer should sign a form. After you have signed the form,
send it by regular mail to

* Johannes Ring
* Simula Research Laboratory
* PO Box 134, 1325 Lysaker
* Norway

Copies of signed consent forms are archived for `authors
<http://fenicsproject.org/pub/copyright/authors>`_ and `institutions
<http://fenicsproject.org/pub/copyright/institutions>`_.

.. _contributing_featured_fenics_stories:

************************************
Contributing Featured FEniCS stories
************************************

The :ref:`Featured_articles` are short descriptions of FEniCS-related
research highlights or new developments.  At any given time, the slide
show on the `FEniCS main page <http://www.fenicsproject.org>`__
highlights a selection of these stories.

The format for a Featured FEniCS story is:

* Graphical highlight -- 366x282 in png/jpg format preferably
* Short summary of research and how FEniCS was used
* Brief information about the author(s)
* Reference list, optionally with a link to repository/Docker image

The text should target a relatively wide audience and should not total
much more than 2-3 paragraphs. Contributions are very welcome --
please send your input to Marie Rognes (meg@simula.no). Plain text
input is perfect, we will format your text to make it web ready and
possibly make minor edits to ensure a consistent style.

.. include linked documents in toctree to avoid Sphinx warning
.. toctree::
    :maxdepth: 1
    :hidden:

    styleguide_cpp
    styleguide_doc
