# Solving the Stokes equations with FEniCS

from fenics import *
from mshr import *

# Define domain
h = 0.25
r = 0.3*h
box = Box(Point(0, 0, 0), Point(1, h, h))
s0 = Sphere(Point(0.3, 0.50*h, 0.50*h), r)
s1 = Sphere(Point(0.5, 0.65*h, 0.65*h), r)
s2 = Sphere(Point(0.7, 0.35*h, 0.35*h), r)
domain = box - s0 - s1 - s2

# Generate mesh
mesh = generate_mesh(domain, 32)
sphere_mesh = generate_mesh(s0 + s1 + s2, 32)

# Define source term
f = Constant((0, 0, 0))

# Define function space
P2 = VectorElement('P', tetrahedron, 2)
P1 = FiniteElement('P', tetrahedron, 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
a = inner(grad(u), grad(v))*dx - p*div(v)*dx + div(u)*q*dx
L = dot(f, v)*dx

# Define boundaries
def inflow(x):
    return near(x[0], 0)

def outflow(x):
    return near(x[0], 1)

def walls(x):
    return near(x[1], 0) or near(x[1], h) or near(x[2], 0) or near(x[2], h)

def spheres(x, on_boundary):
    return on_boundary and not (walls(x) or inflow(x) or outflow(x))

def noslip(x, on_boundary):
    return walls(x) or spheres(x, on_boundary)

# Define boundary conditions
u_D = Expression(('sin(pi*x[1]/h)*sin(pi*x[2]/h)', '0', '0'), h=h, degree=2)
bc0 = DirichletBC(W.sub(0), u_D, inflow)
bc1 = DirichletBC(W.sub(0), (0, 0, 0), noslip)

# Compute solution
w = Function(W)
solve(a == L, w, [bc1, bc0])

# Save solution to file
(u, p) = w.split()
File('u.pvd') << u
File('p.pvd') << p
File('mesh.pvd') << mesh
File('spheres.pvd') << sphere_mesh
