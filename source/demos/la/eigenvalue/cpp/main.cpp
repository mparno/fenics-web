#include <dolfin.h>
#include "Stiffness.h"

using namespace dolfin;

int main()
{
  #ifdef HAS_SLEPC

  // Make sure we use the PETSc backend
  parameters["linear_algebra_backend"] = "PETSc";

  // Create mesh
  Mesh mesh("box_with_dent.xml.gz");

  // Build stiftness matrix
  PETScMatrix A;
  StiffnessMatrix::FunctionSpace V(mesh);
  StiffnessMatrix::BilinearForm a(V, V);
  assemble(a, A);

  // Create eigensolver
  SLEPcEigenSolver eigensolver;

  //Compute all eigenvalues of A x = \lambda x
  eigensolver.solve(A);

  // Extract largest (first) eigenpair
  double r, c;
  PETScVector rx;
  PETScVector cx;
  eigensolver.get_eigenpair(r, c, rx, cx, i);

  std::cout << "Largest eigenvalue: " << r << endl;

  // Initialize function with eigenvector
  Function u(V, rx);

  // Plot eigenfunction
  plot(u);

  }

  #else

    cout << "SLEPc must be installed to run this demo." << endl;

  #endif

  return 0;
}
