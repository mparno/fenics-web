// Copyright (C) 2006-2008 Anders Logg
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
// First added:  2006-02-09
// Last changed: 2013-03-21
//
// This demo solves the Stokes equations, using stabilized
// first order elements for the velocity and pressure. The
// sub domains for the different boundary conditions used
// in this simulation are computed by the demo program in
// src/demo/mesh/subdomains.

#include <dolfin.h>
#include "Stokes.h"

using namespace dolfin;

int main()
{
  // Function for no-slip boundary condition for velocity
  class Noslip : public Expression
  {
  public:

    Noslip() : Expression(2) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      values[0] = 0.0;
      values[1] = 0.0;
    }

  };

  // Function for inflow boundary condition for velocity
  class Inflow : public Expression
  {
  public:

    Inflow() : Expression(2) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      values[0] = -sin(x[1]*DOLFIN_PI);
      values[1] = 0.0;
    }

  };

  // Read mesh and sub domain markers
  Mesh mesh("../dolfin_fine.xml.gz");
  MeshFunction<std::size_t> sub_domains(mesh, "../dolfin_fine_subdomains.xml.gz");

  // Create function space and subspaces
  Stokes::FunctionSpace W(mesh);
  SubSpace W0(W, 0);
  SubSpace W1(W, 1);

  // Create functions for boundary conditions
  Noslip noslip;
  Inflow inflow;
  Constant zero(0.0);

  // No-slip boundary condition for velocity
  DirichletBC bc0(W0, noslip, sub_domains, 0);

  // Inflow boundary condition for velocity
  DirichletBC bc1(W0, inflow, sub_domains, 1);

  // Collect boundary conditions
  std::vector<const DirichletBC*> bcs = {{&bc0, &bc1}};

  // Define variational problem
  Constant f(0.0, 0.0);
  Stokes::BilinearForm a(W, W);
  Stokes::LinearForm L(W);
  L.f = f;
  Function w(W);
  LinearVariationalProblem problem(a, L, w, bcs);

  // Compute solution
  LinearVariationalSolver solver(problem);
  solver.parameters["linear_solver"] = "direct";
  solver.solve();

  // Extract subfunctions
  Function u = w[0];
  Function p = w[1];

  // Save solution in VTK format
  File ufile_pvd("velocity.pvd");
  ufile_pvd << u;
  File pfile_pvd("pressure.pvd");
  pfile_pvd << p;

  // Plot solution
  plot(u);
  plot(p);
  interactive();

  return 0;
}
