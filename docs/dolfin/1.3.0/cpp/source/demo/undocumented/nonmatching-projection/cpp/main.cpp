// Copyright (C) 2009 Garth N. Wells
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
// Modified by Anders Logg, 2011.
//
// First added:  2009-10-10
// Last changed: 2012-11-12
//
// This program demonstrates the L2 projection of a function onto a
// non-matching mesh.

#include <dolfin.h>
#include "P1_projection.h"
#include "P3.h"

using namespace dolfin;

#ifdef HAS_CGAL

class MyExpression : public Expression
{
public:

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(10.0*x[0])*sin(10.0*x[1]);
  }

};

int main()
{
  // Create meshes
  UnitSquareMesh mesh0(16, 16);
  UnitSquareMesh mesh1(64, 64);

  // Create P3 function space
  P3::FunctionSpace V0(mesh0);

  // Interpolate expression into V0
  MyExpression e;
  Function f0(V0);
  f0.interpolate(e);

  // Define variational problem
  P1_projection::FunctionSpace V1(mesh1);
  P1_projection::BilinearForm a(V1, V1);
  P1_projection::LinearForm L(V1, f0);

  // Compute solution
  Function f1(V1);
  solve(a == L, f1);

  // Plot results
  plot(f0);
  plot(f1);
  interactive();

  return 0;
}

#else

int main()
{
  info("DOLFIN must be compiled with CGAL to run this demo.");
  return 0;
}

#endif
