.. _developers_contributing_code:

*****************
Contributing code
*****************

Unless you are a core developer, you will not be able to upload any
code to the main code repositories directly. Instead you will have to
create your own branch and make a merge
request on Bitbucket. If the code is accepted, the patch or branch
will be merged into the main branch by a core developer.

If you are not familiar with `Git <http://git-scm.com/>`__, the
distributed revision control system used for all FEniCS components, a
good starting point is `Bitbucket 101
<https://confluence.atlassian.com/display/BITBUCKET/Bitbucket+101>`__
or the `Git user documentation <http://git-scm.com/documentationl>`__.

.. _before_committing:

Before committing your work
===========================

Before committing any contributions, make sure to test the code
thoroughly. This includes running any unit tests, regression tests
etc. present as part of the code you are modifying. If you are
submitting code for a new feature, it is expected that the new feature
is accompanied by a suitable set of unit tests. You should also review
the :ref:`recommended coding style <developers_writing_code>`.

.. _copyright_consent:

Copyright and licensing consent
===============================

Before your contribution can be accepted into FEniCS, you must sign a
`copyright consent form <http://fenicsproject.org/pub/copyright/forms/>`_.
Ideally, both you and your employer should sign a form. After you have
signed the form, send it by regular mail to

  | Johannes Ring
  | Simula Research Laboratory
  | PO Box 134
  | 1325 Lysaker
  | Norway

Copies of signed consent forms are archived for
`authors <http://fenicsproject.org/pub/copyright/authors>`_
and `institutions <http://fenicsproject.org/pub/copyright/institutions>`_.

Git workflow
============

FEniCS development follows the `gitworkflows
<https://www.kernel.org/pub/software/scm/git/docs/gitworkflows.html>`__
model (with the exception of 'pu' branches). Developers should read up
on the gitworkflows model and understand the role of 'master', 'next'
and 'topic branches'. The same workflow is used by the developers
of `PETSc <http://www.mcs.anl.gov/petsc/>`__. The `PETSc Wiki
<https://bitbucket.org/petsc/petsc/wiki/Home>`__ has some good
information on both Git usage and the gitworkflows model.

A summary of useful Git commands for some common use cases can be found in
the `Git cookbook for FEniCS developers <https://bitbucket.org/fenics-project/dolfin/wiki/Git%20cookbook%20for%20FEniCS%20developers>`__.
