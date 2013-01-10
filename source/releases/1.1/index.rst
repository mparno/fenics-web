.. _release_1_1:

#########################################
Release notes for FEniCS 1.1 (2013-01-07)
#########################################

************
New features
************

The following new features have been introduced in the release
of FEniCS 1.1.0:

* Much improved parallel scaling
* Improve performance of finite element assembly
* Improve performance of mesh topology computation
* Improve performance of mesh building in parallel
* Support for XDMF and HDF5
* Support for restricted function spaces
* Support for for matrix-free linear systems
* Support for GPU-accelerated linear algebra
* Support for mesh generation from constructive solid geometry (CSG) descriptions
* Support for mesh generation from STL
* Support for singular problems with PETSc Krylov solvers
* Change default integer type to std::size_t to handle larger problems
* Make SCOTCH default graph partitioner

For a full list of new features and bug fixes, read the
`ChangeLog <https://launchpad.net/dolfin/1.1.x/1.1.0/>`__

*****************
Interface changes
*****************

We have tried to minimize the number of interface changes between
FEniCS 1.0 and 1.1, but a few changes have been necessary. The key changes are

#. Common mesh classes have been renamed:

   The old classes:

   * Box, Interval, Rectangle, UnitCircle, UnitCube, UnitInterval,
     UnitSphere, UnitSquare, UnitTetrahedron, UnitTriangle

   hav been renamed by adding the suffix 'Mesh'. The new equivalent
   classes are

   * BoxMesh, IntervalMesh, RectangleMesh, UnitCircleMesh,
     UnitCubeMesh, UnitIntervalMesh, UnitSphereMesh, UnitSquareMesh,
     UnitTetrahedronMesh, UnitTriangleMesh

   respectively.

   This change was necessary to differentiate between mesh classes
   the new geometry classes that have been introduced for built-in
   mesh generation.

#. dolfin::uint has been replaced by std::size_t.

   This change was necessary to support large problems with many
   degrees of freedom.

#. The vector of a Function is now returned only as a
   boost::shared_ptr and not as a reference.

#. The AdaptiveLinearVariationalSolvers now take the goal functional
   as a constructor argument (rather than as an argument to solve).
