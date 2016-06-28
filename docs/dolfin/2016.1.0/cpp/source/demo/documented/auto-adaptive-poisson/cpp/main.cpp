// Copyright (C) 2010-2012 Anders Logg and Marie E. Rognes
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
// First added:  2010-08-19
// Last changed: 2012-11-14

#include <dolfin.h>
#include "AdaptivePoisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(5*x[0]);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

int main()
{
  // Create mesh and define function space
  auto mesh = std::make_shared<UnitSquareMesh>(8, 8);
  auto V = std::make_shared<AdaptivePoisson::BilinearForm::TrialSpace>(mesh);

  // Define boundary condition
  auto u0 = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  auto bc = std::make_shared<DirichletBC>(V, u0, boundary);

  // Define variational forms
  auto a = std::make_shared<AdaptivePoisson::BilinearForm>(V, V);
  auto L = std::make_shared<AdaptivePoisson::LinearForm>(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L->f = f;
  L->g = g;

  // Define Function for solution
  auto u = std::make_shared<Function>(V);

  // Define goal functional (quantity of interest)
  auto M = std::make_shared<AdaptivePoisson::GoalFunctional>(mesh);

  // Define error tolerance
  double tol = 1.e-5;

  // Solve equation a = L with respect to u and the given boundary
  // conditions, such that the estimated error (measured in M) is less
  // than tol
  std::vector<std::shared_ptr<const DirichletBC>> bcs({bc});
  auto problem = std::make_shared<LinearVariationalProblem>(a, L, u, bcs);
  AdaptiveLinearVariationalSolver solver(problem, M);
  solver.parameters("error_control")("dual_variational_solver")["linear_solver"]
    = "cg";
  solver.parameters("error_control")("dual_variational_solver")["symmetric"]
    = true;
  solver.solve(tol);

  solver.summary();

  // Plot final solution
  plot(u->root_node(), "Solution on initial mesh");
  plot(u->leaf_node(), "Solution on final mesh");
  interactive();

  return 0;
}
