.. Documentation for the header file dolfin/mesh/Vertex.h

.. _programmers_reference_cpp_mesh_vertex:

Vertex.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Vertex

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
        A Vertex is a MeshEntity of topological dimension 0.

    .. cpp:function:: Point point() const
    
        Return vertex coordinates as a 3D point value

    .. cpp:function:: Vertex(MeshEntity& entity) : MeshEntity(entity.mesh(), 0, entity.index())
    
        Create vertex from mesh entity

    .. cpp:function:: Vertex(const Mesh& mesh, uint index) : MeshEntity(mesh, 0, index)
    
        Create vertex on given mesh

    .. cpp:function:: const double* x() const
    
        Return array of vertex coordinates (const version)

    .. cpp:function:: double x(uint i) const
    
        Return value of vertex coordinate i

    .. cpp:function:: ~Vertex()
    
        Destructor

.. cpp:class:: VertexIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
        A VertexIterator is a MeshEntityIterator of topological dimension 0.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
        A VertexFunction is a MeshFunction of topological dimension 0.

