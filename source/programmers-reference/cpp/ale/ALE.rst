.. Documentation for the header file dolfin/ale/ALE.h

.. _programmers_reference_cpp_ale_ale:

ALE.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: ALE

    This class provides functionality useful for implementation of
    ALE (Arbitrary Lagrangian-Eulerian) methods, in particular
    moving the boundary vertices of a mesh and then interpolating
    the new coordinates for the interior vertices accordingly.

    .. cpp:function:: static void move(Mesh& mesh, BoundaryMesh& new_boundary, dolfin::ALEType method=lagrange)
    
        Move coordinates of mesh according to new boundary coordinates

    .. cpp:function:: static void move(Mesh& mesh, const Function& displacement)
    
        Move coordinates of mesh according to displacement function

    .. cpp:function:: static void move(Mesh& mesh0, Mesh& mesh1, dolfin::ALEType method=lagrange)
    
        Move coordinates of mesh0 according to mesh1 with common global vertices

