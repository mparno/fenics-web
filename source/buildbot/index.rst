.. _fenics_buildbot:

###############
FEniCS buildbot
###############

.. _buildbot_about:

**************
About BuildBot
**************

The `BuildBot <http://www.buildbot.net>`_ is a system to automate the
compile/test cycle required by most software projects to validate code
changes. By automatically rebuilding and testing the tree each time
something has changed, build problems are pinpointed quickly, before
other developers are inconvenienced by the failure. The guilty developer
can be identified and harassed without human intervention. By running
the builds on a variety of platforms, developers who do not have the
facilities to test their changes everywhere before checkin will at least
know shortly afterwards whether they have broken the build or not.

.. _fenics_buildbot_main:

********************
Main FEniCS buildbot
********************

The main FEniCS buildbot builds all :ref:`core components
<about_components_core>` on various platforms in various different ways
each night, runs the test suite, generates :ref:`code coverage reports
<fenics_buildbot_code_coverage>`, and mails any failure to the mailing
lists. To see the current status of the builds, click below.

.. raw:: html
    :file: main_buildbot.inc

.. _fenics_buildbot_personal:

******************
Personal buildbots
******************

We also have a set of personal buildbots that builds DOLFIN from the
personal branches of the core DOLFIN developers. This allows them to
test their code before pushing to the main branch and helps keeping the
main buildbot green.

.. raw:: html
    :file: personal_buildbots.inc

.. _fenics_buildbot_code_coverage:

*************
Code coverage
*************

The buildbot also generates code coverage reports for most of the
:ref:`core components <about_components_core>`. These reports provides
an insight into what parts of the code are executed when the tests are
run.

* `DOLFIN <http://www.fenicsproject.org/coverage/dolfin/lcov/index.html>`__
* `FFC <http://www.fenicsproject.org/coverage/ffc/lcov/index.html>`__
* `FIAT <http://www.fenicsproject.org/coverage/fiat/lcov/index.html>`__
* `Instant <http://www.fenicsproject.org/coverage/instant/lcov/index.html>`__
* `UFL <http://www.fenicsproject.org/coverage/ufl/lcov/index.html>`__
