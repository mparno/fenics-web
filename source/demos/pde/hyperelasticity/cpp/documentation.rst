.. Documentation for the hyperelasticity demo from DOLFIN.

.. _demos_cpp_pde_hyperelasticity:

Hyperelasticity
===============

.. include:: ../common.txt

Implementation
--------------

The implementation is split in two files, a form file containing the
definition of the variational forms expressed in UFL and the solver
which is implemented in a C++ file.

UFL form files
^^^^^^^^^^^^^^

# Function spaces
element = VectorElement("Lagrange", "tetrahedron", 1)

# Test and trial functions
v  = TestFunction(element)      # Test function
du = TrialFunction(element)     # Incremental displacement

# Functions
u  = Coefficient(element)       # Displacement from previous iteration
B  = Coefficient(element)       # Body force per unit mass
T  = Coefficient(element)       # Traction force on the boundary

# Kinematics
I = Identity(element.cell().d)  # Identity tensor
F = I + grad(u)                 # Deformation gradient
C = F.T*F                       # Right Cauchy-Green tensor

# Invariants of deformation tensors
Ic = tr(C)
J  = det(F)

# Elasticity parameters
mu    = Constant("tetrahedron")
lmbda = Constant("tetrahedron")

# Stored strain energy density (compressible neo-Hookean model)
psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Total potential energy
Pi = psi*dx - inner(B, u)*dx - inner(T, u)*ds

# First variation of Pi (directional derivative about u in the direction of v)
L = derivative(Pi, u, v)

# Compute Jacobian of L
a = derivative(L, u, du)

C++ code
^^^^^^^^

.. code-block:: cpp

  #include <dolfin.h>
  #include "HyperElasticity.h"

.. code-block:: cpp

  using namespace dolfin;

.. code-block:: cpp

  // Sub domain for clamp at left end
  class Left : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return (std::abs(x[0]) < DOLFIN_EPS) && on_boundary;
    }
  };

.. code-block:: cpp

  // Sub domain for rotation at right end
  class Right : public SubDomain
  {
    bool inside(const Array<double>& x, bool on_boundary) const
    {
      return (std::abs(x[0] - 1.0) < DOLFIN_EPS) && on_boundary;
    }
  };

.. code-block:: cpp

  // Dirichlet boundary condition for clamp at left end
  class Clamp : public Expression
  {
  public:

    Clamp() : Expression(3) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
    }

  };

.. code-block:: cpp

  // Dirichlet boundary condition for rotation at right end
  class Rotation : public Expression
  {
  public:

    Rotation() : Expression(3) {}

    void eval(Array<double>& values, const Array<double>& x) const
    {
      const double scale = 0.5;

      // Center of rotation
      const double y0 = 0.5;
      const double z0 = 0.5;

      // Large angle of rotation (60 degrees)
      double theta = 1.04719755;

      // New coordinates
      double y = y0 + (x[1] - y0)*cos(theta) - (x[2] - z0)*sin(theta);
      double z = z0 + (x[1] - y0)*sin(theta) + (x[2] - z0)*cos(theta);

      // Rotate at right end
      values[0] = 0.0;
      values[1] = scale*(y - x[1]);
      values[2] = scale*(z - x[2]);
    }
  };

.. code-block:: cpp

  int main()
  {

.. code-block:: cpp

  // Create mesh and define function space
  UnitCube mesh (16, 16, 16);
  HyperElasticity::FunctionSpace V(mesh);

.. code-block:: cpp

  // Define Dirichlet boundaries
  Left left;
  Right right;

  // Define Dirichlet boundary functions
  Clamp c;
  Rotation r;

  // Create Dirichlet boundary conditions
  DirichletBC bcl(V, c, left);
  DirichletBC bcr(V, r, right);
  std::vector<const BoundaryCondition*> bcs;
  bcs.push_back(&bcl); bcs.push_back(&bcr);

.. code-block:: cpp

  // Define source and boundary traction functions
  Constant B(0.0, -0.5, 0.0);
  Constant T(0.1,  0.0, 0.0);

.. code-block:: cpp

  // Define solution function
  Function u(V);

.. code-block:: cpp

  // Set material parameters
  const double E  = 10.0;
  const double nu = 0.3;
  Constant mu(E/(2*(1 + nu)));
  Constant lambda(E*nu/((1 + nu)*(1 - 2*nu)));

.. code-block:: cpp

  // Create forms
  HyperElasticity::BilinearForm a(V, V);
  a.mu = mu; a.lmbda = lambda; a.u = u;
  HyperElasticity::LinearForm L(V);
  L.mu = mu; L.lmbda = lambda; L.B = B; L.T = T; L.u = u;

.. code-block:: cpp

  // Solve nonlinear variational problem
  VariationalProblem problem(a, L, bcs, true);
  problem.solve(u);

.. code-block:: cpp

  // Save solution in VTK format
  File file("displacement.pvd");
  file << u;

  // Plot solution
  plot(u);

  return 0;
}

Complete code
-------------

.. literalinclude:: main.cpp
   :start-after: // Begin demo
   :language: c++

