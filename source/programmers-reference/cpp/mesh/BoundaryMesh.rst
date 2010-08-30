.. Documentation for the header file dolfin/mesh/BoundaryMesh.h

.. _programmers_reference_cpp_mesh_boundarymesh:

BoundaryMesh.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: BoundaryMesh

    *Parent class*
    
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

    .. cpp:function:: ~BoundaryMesh()
    
        Destructor

