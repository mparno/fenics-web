.. _release_1_4:

#########################################
Release notes for FEniCS 1.4 (2014-06-02)
#########################################

************
New features
************

The following new features have been introduced in the release
of FEniCS 1.4.0:

* Merge UFC into FFC; UFC is now distributed as part of FFC
* Add support for VTK 6
* Reduce memory usage and increase speed of mesh topology computation
* Speed up of dof reordering for mixed spaces including global dofs
* Add new `fenics` module, as an alias for `dolfin`
* Remove Boost.MPI dependency
* Add experimental support for new integral type `custom_integral`
* Improve svg rendering of cells and sobolev spaces in iPython notebook
* Improve restriction handling; restricting continuous coefficients and constants is now optional
* Introduce several new geometric quantities in the form language
* Improve notation for domains in the form language
* Fix several bugs related to facet integrals over manifolds

In addition to the changes listed above, a large number of issues
have been fixed since the release of FEniCS 1.3.0.

For a full list of new features and bug fixes, read the
ChangeLogs for
`DOLFIN <http://fenicsproject.org/pub/software/dolfin/ChangeLog>`__,
`FFC <http://fenicsproject.org/pub/software/ffc/ChangeLog>`__, and
`UFL <http://fenicsproject.org/pub/software/ufl/ChangeLog>`__.

*****************
Interface changes
*****************

The following changes have been made to the FEniCS interface between
versions 1.3 and 1.4:

* The mesh classes `UnitSquare` and friends have now been removed (replaced by `UnitSquareMesh` etc)
* `size_t` should be used throughout in place of `uint`
