// Copyright (C) 2007-2011 Kristian B. Oelgaard, Anders Logg and Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2007-06-29
// Last changed: 2013-03-21
//
// Steady state advection-diffusion equation, discontinuous
// formulation using full upwinding.

#include <dolfin.h>

#include "AdvectionDiffusion.h"
#include "Projection.h"
#include "Velocity.h"

using namespace dolfin;

// Dirichlet boundary condition
class BC : public Expression
{
public:

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(DOLFIN_PI*5.0*x[1]);
  }

};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return std::abs(x[0] - 1.0) < DOLFIN_EPS && on_boundary;
  }
};

int main(int argc, char *argv[])
{
  // FIXME: Make mesh ghosted
  parameters["ghost_mode"] = "shared_facet";

  // Read simple velocity field (-1.0, -0.4) defined on a 64x64 unit square
  // mesh and a quadratic vector Lagrange element

  // Read mesh
  Mesh mesh("../unitsquare_64_64.xml.gz");

  // Create velocity FunctionSpace
  Velocity::FunctionSpace V_u(mesh);

  // Create velocity function
  Function u(V_u, "../unitsquare_64_64_velocity.xml.gz");

  // Diffusivity
  Constant c(0.0);

  //Source term
  Constant f(0.0);

  // Penalty parameter
  Constant alpha(5.0);

  // Create function space
  AdvectionDiffusion::FunctionSpace V(mesh);

  // Create forms and attach functions
  AdvectionDiffusion::BilinearForm a(V, V);
  a.u = u; a.kappa = c; a.alpha = alpha;
  AdvectionDiffusion::LinearForm L(V);
  L.f = f;

  // Set up boundary condition (apply strong BCs)
  BC g;
  DirichletBoundary boundary;
  DirichletBC bc(V, g, boundary, "geometric");

  // Solution function
  Function phi_h(V);

  // Assemble and apply boundary conditions
  Matrix A;
  Vector b;
  assemble(A, a);
  assemble(b, L);
  bc.apply(A, b);

  // Solve system
  solve(A, *phi_h.vector(), b);

  // Define variational problem
  Projection::FunctionSpace Vp(mesh);
  Projection::BilinearForm ap(Vp, Vp);
  Projection::LinearForm Lp(Vp);
  Lp.phi0 = phi_h;

  // Compute solution
  Function phi_p(Vp);
  solve(ap == Lp, phi_p);

  // Save projected solution in VTK format
  File file("temperature.pvd");
  file << phi_p;

  // Plot projected solution
  plot(phi_h);
  interactive();
}
