
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


    .. cpp:function:: BoundaryMesh()
    
        Create an empty boundary mesh


    .. cpp:function:: BoundaryMesh(const Mesh& mesh, bool order=true)
    
        Create boundary mesh from given mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Another :cpp:class:`Mesh` object.
            order (bool)
                Optional argument which can be used to control whether
                or not the boundary mesh should be ordered according
                to the UFC ordering convention. If set to false, the
                boundary mesh will be ordered with right-oriented
                facets (outward-pointing unit normals). The default
                value is true.


    .. cpp:function:: void init_exterior_boundary(const Mesh& mesh)
    
        Initialize exterior boundary of given mesh


    .. cpp:function:: void init_interior_boundary(const Mesh& mesh)
    
        Initialize interior boundary of given mesh


    .. cpp:function:: const MeshFunction<std::size_t>& cell_map() const
    
        Get cell mapping from the boundary mesh to the original full mesh


    .. cpp:function:: MeshFunction<std::size_t>& vertex_map()
    
        Get vertex mapping from the boundary mesh to the original full mesh


