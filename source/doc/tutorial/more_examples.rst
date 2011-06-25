.. Automatically generated reST file from Doconce source
   (http://code.google.com/p/doconce/)


More Examples
=============

Many more topics could be treated in a FEniCS tutorial, e.g., how
to solve systems of PDEs, how to work with mixed finite element
methods, how to create more complicated meshes and mark boundaries,
and how to create more advanced visualizations.  However, to limit the
size of this tutorial, the examples end here.  There are, fortunately,
a rich set of examples coming with the DOLFIN source code. Go to
``dolfin/demo``. The subdirectory ``pde`` contains many examples
on solving PDEs:

  * the advection-diffusion equation (``advection-diffusion``),

  * the Cahn-Hilliard equation (``cahn-hilliard``),

  * the equation of linear elasticity (``elasticity``)

and hyperelasticity (``hyperelasticity``),
  * the Poisson equation with a variable tensor coefficient

(``tensor-weighted-poisson``),
  * mixed finite elements for the Poisson

equation (``mixed-poisson``),
  * the Stokes problem of fluid flow (``stokes``),

  * an eigenvalue problem arising from electromagnetic

waveguide problem with N\'{e}d\'{e}lec elements.

Moreover, the ``dg`` subdirectory contains demonstrations of
applying discontinuous Galerkin methods to
the advection-diffusion, Poisson, and
Biharmonic equations.
There also exists an example on how to compute
functionals over subsets of the mesh (``lift-drag``).

The ``demo/mesh`` directory contains examples on moving a mesh
(``ale``), computing intersections (``intersection``),
mesh refinement (``refinement``), and creating separate subdomain
meshes from a common parent mesh (``submesh``).

The ``cbc.solve`` suite of applications is under development and
will contain Navier--Stokes solvers and large-strain elasticity
solvers.  The ``cbc.rans`` suite will in particular contain several
Navier--Stokes solvers in combination with a range of PDEs arising in
various turbulence models.
