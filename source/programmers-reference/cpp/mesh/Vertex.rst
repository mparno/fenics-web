.. Documentation for the header file dolfin/mesh/Vertex.h

.. _programmers_reference_cpp_mesh_vertex:

Vertex.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Vertex

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
    A Vertex is a MeshEntity of topological dimension 0.

    .. cpp:function:: Vertex(const Mesh& mesh, uint index)
    
        Create vertex on given mesh

    .. cpp:function:: Vertex(MeshEntity& entity)
    
        Create vertex from mesh entity

    .. cpp:function:: double x(uint i) const
    
        Return value of vertex coordinate i

    .. cpp:function:: Point point() const
    
        Return vertex coordinates as a 3D point value

    .. cpp:function:: const double* x() const
    
        Return array of vertex coordinates (const version)

.. cpp:class:: VertexIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
    A VertexIterator is a MeshEntityIterator of topological dimension 0.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
    A VertexFunction is a MeshFunction of topological dimension 0.

