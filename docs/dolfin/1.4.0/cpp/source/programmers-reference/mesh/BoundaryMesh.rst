
.. Documentation for the header file dolfin/mesh/BoundaryMesh.h

.. _programmers_reference_cpp_mesh_boundarymesh:

BoundaryMesh.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BoundaryMesh

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    A BoundaryMesh is a mesh over the boundary of some given mesh.
    The cells of the boundary mesh (facets of the original mesh) are
    oriented to produce outward pointing normals relative to the
    original mesh.


    .. cpp:function:: BoundaryMesh(const Mesh& mesh, std::string type, bool order=true)
    
        Create boundary mesh from given mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Another :cpp:class:`Mesh` object.
            type (_std::string_)
                The type of BoundaryMesh, which can be "exterior",
                "interior" or "local". "exterior" is the globally
                external boundary, "interior" is the inter-process mesh
                and "local" is the boudary of the local (this process)
                mesh.
            order (bool)
                Optional argument which can be used to control whether
                or not the boundary mesh should be ordered according
                to the UFC ordering convention. If set to false, the
                boundary mesh will be ordered with right-oriented
                facets (outward-pointing unit normals). The default
                value is true.


    .. cpp:function:: MeshFunction<std::size_t>& entity_map(std::size_t d)
    
        Get index map for entities of dimension d in the boundary mesh
        to the entity in the original full mesh


    .. cpp:function:: const MeshFunction<std::size_t>& entity_map(std::size_t d) const
    
        Get index map for entities of dimension d in the boundary mesh
        to the entity in the original full mesh (const version)


