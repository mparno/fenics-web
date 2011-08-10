
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


    .. cpp:function:: BoundaryMesh()
    
        Create an empty boundary mesh


    .. cpp:function:: BoundaryMesh(const Mesh& mesh)
    
        Create (interior) boundary mesh from given mesh


    .. cpp:function:: void init_exterior_boundary(const Mesh& mesh)
    
        Initialize exterior boundary of given mesh


    .. cpp:function:: void init_interior_boundary(const Mesh& mesh)
    
        Initialize interior boundary of given mesh


    .. cpp:function:: const MeshFunction<unsigned int>& cell_map() const
    
        Get cell mapping from the boundary mesh to the original full mesh


    .. cpp:function:: MeshFunction<unsigned int>& vertex_map()
    
        Get vertex mapping from the boundary mesh to the original full mesh


