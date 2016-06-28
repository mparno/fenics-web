
.. Documentation for the header file dolfin/mesh/Vertex.h

.. _programmers_reference_cpp_mesh_vertex:

Vertex.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Vertex

    *Parent class(es)*
    
        * :cpp:class:`MeshEntity`
        
    A Vertex is a MeshEntity of topological dimension 0.


    .. cpp:function:: Vertex(const Mesh& mesh, std::size_t index)
    
        Create vertex on given mesh


    .. cpp:function:: Vertex(MeshEntity& entity)
    
        Create vertex from mesh entity


    .. cpp:function:: double x(std::size_t i) const
    
        Return value of vertex coordinate i


    .. cpp:function:: Point point() const
    
        Return vertex coordinates as a 3D point value


    .. cpp:function:: const double* x() const
    
        Return array of vertex coordinates (const version)


.. cpp:class:: VertexFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A VertexFunction is a MeshFunction of topological dimension 0.


