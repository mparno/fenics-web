.. General introduction to the FEniCS documentation effort.
   This is where we explain the main idea and structure of the docs.

.. _introduction:

############
Introduction
############

This is the documentation for the
`FEniCS project <http://fenics.org/wiki/FEniCS_Project>`_.
The FEniCS project consists of many
`software components <http://fenics.org/wiki/Projects>`_, however, the
focus of this documentation is on `DOLFIN <http://fenics.org/wiki/DOLFIN>`_
which is the main interface to FEniCS.
The documentation comprises three main parts:

* :ref:`tutorial_index`, which focuses on getting users started solving
  partial differential equations (PDEs) right away.
  The tutorial is self-contained and shows how to use FEniCS tools to setup
  and solve problems.
  Although the tutorial starts with a simple problem, the complexity is soon
  increased to demonstrate how intuitive FEniCS software can be put to work.
  This is the ideal place to start for newcomers to FEniCS and people who are
  simply curious to learn what exactly FEniCS is all about.

* :ref:`demos_index`, which is a collection of *ready to run code examples*
  that demonstrate a subset of FEniCS features.
  This typically involves solving a particular PDE but can just as well be
  demonstrating functionality pertaining to mesh manipulation, linear algebra
  etc.
  Another purpose of the demos is to serve as code templates for FEniCS users
  such that solvers does not have to be written from scratch since it is
  usually possible to find a demo which bears resemblance to the problem at
  hand and modify it accordingly.

* :ref:`programmers_reference_index`, which contains detailed documentation
  of the functions and classes in DOLFIN.


All three parts come in a ``C++`` and a ``Python`` version and they can be
accessed online as ``HTML`` or downloaded as ``PDF`` files.

Any suggestions for improving FEniCS software or the documentation is more
than welcome, please see :ref:`contributing` on how to do this.

