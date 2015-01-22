.. _release_1_5:

#########################################
Release notes for FEniCS 1.5 (2015-01-12)
#########################################

************
New features
************

The following new features have been introduced in the release
of FEniCS 1.5.0:

* Parallel assembly of DG methods
* Parallel adaptive mesh refinement in both 2D and 3D
* Ghosted meshes
* Full support for linear algebra backends with 64-bit integers
* Significant memory reduction in dofmap storage
* Switch to local (process-wise) indexing for dof indices
* Experimental support for assembly over multiple non-matching meshes (multimesh)
* Support for assembly of `custom_integral`
* Move mesh generation functionality to `mshr`
* Fixes for petsc4py and slepc4py
* Remove CGAL dependency
* CMake 3 compatibility
* PETSc 3.5 compatibility, require version >= 3.3
* Python 3 experimental support, require version >= 2.7

In addition to the changes listed above, a large number of issues
have been fixed since the release of FEniCS 1.4.0.

For a full list of new features and bug fixes, read the
ChangeLogs for
`DOLFIN <http://fenicsproject.org/pub/software/dolfin/ChangeLog>`__,
`FFC <http://fenicsproject.org/pub/software/ffc/ChangeLog>`__, and
`UFL <http://fenicsproject.org/pub/software/ufl/ChangeLog>`__.

*****************
Interface changes
*****************

The following changes have been made to the FEniCS interface between
versions 1.4 and 1.5:

* The classes `FacetNormal` and friends now require a mesh instead of a cell.
* The signatures for assemble() have been simplified to not require/accept `coefficients`, `cells`, and `ccommon_cell`.
