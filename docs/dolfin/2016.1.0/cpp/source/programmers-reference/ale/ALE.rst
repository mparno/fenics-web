
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


    .. cpp:function:: static std::shared_ptr<MeshDisplacement> move(std::shared_ptr<Mesh> mesh, const BoundaryMesh& new_boundary)
    
        Move coordinates of mesh according to new boundary coordinates.
        Works only for affine meshes.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The affine mesh to move.
            boundary (:cpp:class:`BoundaryMesh`)
                An affine mesh containing just the boundary cells.
        
        *Returns*
            MeshDisplacement
                Displacement encapsulated in Expression subclass
                MeshDisplacement.


    .. cpp:function:: static std::shared_ptr<MeshDisplacement> move(std::shared_ptr<Mesh> mesh0, const Mesh& mesh1)
    
        *Arguments*
            mesh0 (:cpp:class:`Mesh`)
                The affine mesh to move.
            mesh1 (:cpp:class:`Mesh`)
                The affine mesh to be fit.
        
        *Returns*
            MeshDisplacement
                Displacement encapsulated in Expression subclass


    .. cpp:function:: static void move(Mesh& mesh, const GenericFunction& displacement)
    
        Move coordinates of mesh according to displacement function.
        This works only for affine meshes.
        
        NOTE: This cannot be implemented for higher-order geometries
              as there is no way of constructing function space for
              position unless supplied as an argument.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The affine mesh to move.
            displacement (:cpp:class:`GenericFunction`)
                A vectorial generic function.


    .. cpp:function:: static void move(Mesh& mesh, const Function& displacement)
    
        Move coordinates of mesh according to displacement function.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to move.
            displacement (:cpp:class:`Function`)
                A vectorial Lagrange function of matching degree.


