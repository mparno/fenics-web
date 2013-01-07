.. _release_1_1:

############################
Release notes for FEniCS 1.1
############################

We have tried to minimize the number of interface changes between
FEniCS 1.0 and 1.1, but a few has been necessary. The key changes are

#. Common mesh classes have been renamed:

   The old classes:

   * Box, Interval, Rectangle, UnitCircle, UnitCube, UnitInterval,
     UnitSphere, UnitSquare, UnitTetrahedron, UnitTriangle

   has been renamed by adding the suffix 'Mesh'. The new equivalent
   classes are

   * BoxMesh, IntervalMesh, RectangleMesh, UnitCircleMesh,
     UnitCubeMesh, UnitIntervalMesh, UnitSphereMesh, UnitSquareMesh,
     UnitTetrahedronMesh, UnitTriangleMesh

   respectively.

   Note in particular that the classes Box and Rectangle now exists as
   csg shape classes. For the rest of the classes, deprecation classes
   with the old are still present.

#. dolfin::uint has been removed and size_t has been introduced

#. The vector of a Function is now returned only as a
   boost::shared_ptr and not as a reference.

#. The AdaptiveLinearVariationalSolvers now take the goal functional
   as a constructor argument (rather than as an argument to solve).
