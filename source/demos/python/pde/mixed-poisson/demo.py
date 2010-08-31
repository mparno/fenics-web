"""This demo program solves the mixed formulation of Poisson's
equation:

    sigma + grad(u) = 0
         div(sigma) = f

The corresponding weak (variational problem)

    <tau, sigma> - <div(tau), u> = 0       for all tau
                 <v, div(sigma)> = <w, f>  for all v

is solved using BDM (Brezzi-Douglas-Marini) elements of degree q (tau,
sigma) and DG (discontinuous Galerkin) elements of degree q - 1 for
(w, u).

Original implementation: ../cpp/main.cpp by Anders Logg and Marie Rognes
"""

__author__    = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__      = "2007-11-14 -- 2008-12-19"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__   = "GNU LGPL Version 2.1"

# Modified by Marie E. Rognes 2010

# Begin demo

from dolfin import *

# Create mesh
mesh = UnitSquare(32, 32)

# Define function spaces and mixed (product) space
BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)
V = BDM * DG

# Define trial and test functions
(sigma, u) = TrialFunctions(V)
(tau, v) = TestFunctions(V)

# Define source function
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")

# Define variational form
a = (dot(sigma, tau) - u*div(tau) + div(sigma)*v)*dx
L = v*f*dx

# Define essential boundary (y = 0 or y = 1)
def on_boundary(x):
    return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

# Define function G such that G \cdot n = g
class Flux(Expression):
    def eval_data(self, values, data):
        g = - sin(5*data.x()[0])
        values[0] = g*data.normal()[0]
        values[1] = g*data.normal()[1]
    def rank(self):
        return 1
    def dim(self):
        return 2
G = Flux()
bc = DirichletBC(V.sub(0), G, on_boundary)

# Compute solution and plot u
problem = VariationalProblem(a, L, bc)
(sigma, u) = problem.solve().split()
plot(u, interactive=True)
