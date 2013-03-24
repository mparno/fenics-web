
.. Documentation for the header file dolfin/ale/ALE.h

.. _programmers_reference_cpp_ale_ale:

ALE.h
=====

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ALE

    This class provides functionality useful for implementation of
    ALE (Arbitrary Lagrangian-Eulerian) methods, in particular
    moving the boundary vertices of a mesh and then interpolating
    the new coordinates for the interior vertices accordingly.


    .. cpp:function:: static boost::shared_ptr<MeshDisplacement> move(Mesh& mesh, const BoundaryMesh& new_boundary)
    
        Move coordinates of mesh according to new boundary coordinates.
        Returns displacement (encapsulated in Expression subclass MeshDisplacement)


    .. cpp:function:: static boost::shared_ptr<MeshDisplacement> move(Mesh& mesh0, const Mesh& mesh1)
    
        Move coordinates of mesh0 according to mesh1 with common global vertices.
        Returns displacement (encapsulated in Expression subclass MeshDisplacement)


    .. cpp:function:: static void move(Mesh& mesh, const GenericFunction& displacement)
    
        Move coordinates of mesh according to displacement function


