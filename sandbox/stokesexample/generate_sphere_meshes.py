# Generate sphere meshes to support the visualization.
# The standard FEniCS mesh will not be oriented which
# sqrews up the visualization in Paraview. This script
# is a hack so we can generate a pretty picture.
#
# Orienting the meshes actually doesn't seem to help.
# Paraview seems to color the faces based on their area
# or some other property unrelated to their orientation.
# Generating very fine meshes for the spheres instead
# which at least makes them look smooth.

from fenics import *
from mshr import *
import numpy

def orient_mesh(mesh):

    # Find center of mesh
    x = mesh.coordinates()
    center = Point(numpy.mean(x[:, 0]), numpy.mean(x[:, 1]), numpy.mean(x[:, 2]))

    # Open new mesh for editing
    editor = MeshEditor()
    new_mesh = Mesh()
    editor.open(new_mesh, 2, 3)
    editor.init_vertices(mesh.num_vertices())
    editor.init_cells(mesh.num_cells())

    # Add all vertices
    for i, v in enumerate(vertices(mesh)):
        editor.add_vertex(i, v.point())

    # Add all cells and fix orientation
    g = mesh.geometry()
    for i, c in enumerate(cells(mesh)):
        v = [v for v in vertices(c)]
        w1 = v[1].point() - v[0].point()
        w2 = v[2].point() - v[0].point()
        orientation = w1.cross(w2).dot(v[0].point() - center)
        if orientation > 0:
           print 'orientation = 1'
           editor.add_cell(i, v[0].index(), v[1].index(), v[2].index())
        else:
           print 'orientation = -1'
           editor.add_cell(i, v[0].index(), v[2].index(), v[1].index())

    # Close mesh editor
    editor.close()

    return new_mesh

# Define domain
h = 0.25
r = 0.3*h
box = Box(Point(0, 0, 0), Point(1, h, h))
s0 = Sphere(Point(0.3, 0.50*h, 0.50*h), r)
s1 = Sphere(Point(0.5, 0.65*h, 0.65*h), r)
s2 = Sphere(Point(0.7, 0.35*h, 0.35*h), r)

# Generate meshes
n = 40
s0 = BoundaryMesh(generate_mesh(s0, n), 'exterior')
s1 = BoundaryMesh(generate_mesh(s1, n), 'exterior')
s2 = BoundaryMesh(generate_mesh(s2, n), 'exterior')

# Orient meshes
s0 = orient_mesh(s0)
s1 = orient_mesh(s1)
s2 = orient_mesh(s2)

# Save meshes to file
File('s0.pvd') << s0
File('s1.pvd') << s1
File('s2.pvd') << s2
