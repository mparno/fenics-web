"""
This program illustrates basic use of the SLEPc eigenvalue solver for
a standard eigenvalue problem."""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-11-28 -- 2009-10-09"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU LGPL Version 2.1"

# Modified by Anders Logg, 2008.
# Modified by Marie Rognes, 2009.

# Begin demo

from dolfin import *

# Test for PETSc and SLEPc
if not has_la_backend("PETSc"):
    print "DOLFIN has not been configured with PETSc. Exiting."
    exit()

if not has_slepc():
    print "DOLFIN has not been configured with SLEPc. Exiting."
    exit()

# Make sure we use the PETSc backend
parameters["linear_algebra_backend"] = "PETSc"

# Define mesh, function space and basis functions
mesh = UnitSquare(4, 4)
V = FunctionSpace(mesh, "CG", 1)
v = TestFunction(V)
u = TrialFunction(V)

# Define form for stiffness matrix
a = dot(grad(v), grad(u))*dx

# Assemble stiffness form
A = PETScMatrix()
assemble(a, tensor=A)

# Compute all eigenvalues of A x = \lambda x
esolver = SLEPcEigenSolver()
esolver.solve(A)

# Get largest (first) eigenpair
a, b, x, y = esolver.get_eigenpair(0)
print "Largest eigenvalue: ", a

# Initialize function with eigenvector
u = Function(V, x)

# Plot eigenfunction
plot(u)
interactive()
