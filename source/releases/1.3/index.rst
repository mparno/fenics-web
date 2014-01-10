.. _release_1_3:

#########################################
Release notes for FEniCS 1.3 (2014-01-07)
#########################################

************
New features
************

The following new features have been introduced in the release
of FEniCS 1.3.0:

* Assignment of sub functions
* Mesh distance queries in 3D
* Solving of local ODEs at vertices
* Runge-Kutta and multistage ODE solvers
* New operators ``cell_avg`` and ``facet_avg``
* New built-in library for computational geometry (bounding box trees)
* Updates for NumPy 1.7

For a full list of new features and bug fixes, read the
`ChangeLog <http://fenicsproject.org/pub/software/dolfin/ChangeLog>`__.

*****************
Interface changes
*****************

The following changes have been made to the FEniCS interface between
versions 1.2 and 1.3:

* The function ``compiled_subdomains`` has been replaced by the class ``CompiledSubDomain``.
