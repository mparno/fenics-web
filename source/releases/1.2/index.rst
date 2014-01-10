.. _release_1_2:

#########################################
Release notes for FEniCS 1.2 (2013-03-24)
#########################################

************
New features
************

The following new features have been introduced in the release
of FEniCS 1.2.0:

* Support for solving PDEs on manifolds
* New and improved implementation of periodic boundary conditions

For a full list of new features and bug fixes, read the
`ChangeLog <http://fenicsproject.org/pub/software/dolfin/ChangeLog>`__.

*****************
Interface changes
*****************

The following changes have been made to the FEniCS interface between
versions 1.1 and 1.2:

* Optional arguments to assemblers have been removed
* The class SymmetricAssembler has been removed
* Domains for assemblers can now only be attached to forms
