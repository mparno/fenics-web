from fenics import *

# Import the mesh
mesh = Mesh("klein.xml.gz")

# Define the function space of choice as usual
V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

# Define the variational form as usual, but pay attention to the
# definition of derivatives and measures
f = Expression("x[0]*x[1]*x[2]", degree=3)
c = Constant(0.03)
a = c*u*v*dx + inner(grad(u), grad(v))*dx
L = f*v*dx

# Solve and plot as usual
u = Function(V)
solve(a == L, u)

plot(mesh, title="Gray's Klein bottle")
plot(u, title="u", interactive=True)
