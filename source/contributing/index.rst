.. Developer information.

.. _developers:

############
Contributing
############

FEniCS is a free/open source project and everyone is invited to
contribute. This page contains information for prospective FEniCS
developers, including an overview of the organization of the FEniCS
project, how to write code and documentation, and how to use tools
like Bitbucket and Git.

************
Organization
************

FEniCS is organized as a collection of interoperable components that
together form the FEniCS Project. Each component is developed by one
or more authors. This means that each component can be developed at
its own pace. At the same time, we strive to make periodic and
coordinated releases of all components to ensure interoperability
between the components.

The FEniCS Project uses `Bitbucket <https://bitbucket.org>`_ as its
main development platform. All FEniCS projects are collected under
the `FEniCS Project on Bitbucket
<https://bitbucket.org/fenics-project>`_.

*************************
Obtaining the source code
*************************

FEniCS uses `Bitbucket <http://bitbucket.org>`__ for hosting
code. Each FEniCS component has a `Git <http://git-scm.com/>`_
repository on Bitbucket that contains all source code (including the
entire development history). The repositories are readable for
everyone, but write access is only granted to the members of the core
teams.

======================================
Accessing the development repositories
======================================

To access the development repositories, you first need to install the
revision control system Git. Visit the `Git web page
<http://git-scm.com/>`__ for instructions on how to install Git on
your operating system.

Once Git has been installed, you can access the development repository
of any FEniCS project by the ``git`` command. For example, to check
out the source code for DOLFIN, simply issue the following command::

    git clone https://bitbucket.org/fenics-project/dolfin.git

or if you have an ssh key associated with your Bitbucket account::

    git clone git@bitbucket.org:fenics-project/dolfin.git

========================
Notifications of updates
========================

Developers should subscribe to notifications of changes made to the
source code by visiting the repository on Bitbucket and clicking the
'follow' button.

Links to the source repositories for FEniCS projects can be found on
the `FEniCS Bitbucket page <https://bitbucket.org/fenics-project>`__.

======================
FEniCS Developer Tools
======================

Developers should take a look at the `FEniCS Developer Tools
<https://bitbucket.org/fenics-project/fenics-developer-tools>`__
repository. This contains scripts that are very useful for developers
in setting up a good development environment for FEniCS. In
particular, consider using the scripts ``fenics-install-all.sh`` and
``fenics-install-component.sh``.

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

.. _before_committing:

===========================
Before committing your work
===========================

Before committing any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying. If you are
submitting code for a new feature, it is expected that the new feature
is accompanied by a suitable set of unit tests. You should also review
the :ref:`recommended coding style <developers_writing_code>`.

.. _copyright_consent:

===============================
Copyright and licensing consent
===============================

Before your contribution can be accepted into FEniCS, you must sign a
`copyright consent form
<http://fenicsproject.org/pub/copyright/forms/>`_.  Ideally, both you
and your employer should sign a form. After you have signed the form,
send it by regular mail to

  | Johannes Ring
  | Simula Research Laboratory
  | PO Box 134
  | 1325 Lysaker
  | Norway

Copies of signed consent forms are archived for `authors
<http://fenicsproject.org/pub/copyright/authors>`_ and `institutions
<http://fenicsproject.org/pub/copyright/institutions>`_.

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

*************************************
Contributing FEniCS Featured articles
*************************************

The :ref:`Featured_articles` are short descriptions of FEniCS
highlights: new FEniCS features, new applications of FEniCS or current
FEniCS Events. At any given time, the slide show on the `FEniCS main
page <http://www.fenicsproject.org>`__. highlights a selection of
these articles.

All FEniCS developers and users are encouraged to contribute Featured
articles. Detailed instructions are given here: `How to contribute a
FEniCS Featured article
<https://bitbucket.org/snippets/fenics-project/LdRGq/how-to-create-a-new-fenics-featured>`__.

.. _developers_writing_code:

************
Writing code
************

We ask all developers and contributors to adhere to a common style
guide. This makes the job easier for maintainers who need to review,
edit and maintain the FEniCS code base.

Patches that don't follow the correct coding style will likely be
rejected, as the maintainer responsible for reviewing the patch must
otherwise make additional efforts to edit the patch to follow the
coding style.

===============
C++ style guide
===============

Guidelines for writing C++ code are given in the
:ref:`developers_styleguide_cpp`.

==================
Python style guide
==================

The FEniCS coding style for Python code adheres to the `PEP-8 style
guide <http://www.python.org/dev/peps/pep-0008/>`_ although it is not
strictly enforced.

=========================
Documentation style guide
=========================

The FEniCS documentation is generated by `Sphinx
<http://sphinx.pocoo.org/index.html>`_ and uses `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_ (reST) as the markup
language.  Good starting points are `reStructuredText primer
<http://sphinx.pocoo.org/rest.html>`_ and `Sphinx Markup Constructs
<http://sphinx.pocoo.org/markup/index.html>`_.

Further guidelines are given in the :ref:`styleguide_documentation`.

.. include linked documents in toctree to avoid Sphinx warning
.. toctree::
    :maxdepth: 1
    :hidden:

    styleguide_cpp
    styleguide_doc
