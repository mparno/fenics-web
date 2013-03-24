
.. Documentation for the header file dolfin/mesh/Edge.h

.. _programmers_reference_cpp_mesh_edge:

Edge.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Edge

    *Parent class(es)*
    
        * :cpp:class:`MeshEntity`
        
    An Edge is a :cpp:class:`MeshEntity` of topological dimension 1.


    .. cpp:function:: Edge(const Mesh& mesh, std::size_t index)
    
        Create edge on given mesh
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            index (std::size_t)
                Index of the edge.


    .. cpp:function:: Edge(MeshEntity& entity)
    
        Create edge from mesh entity
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The mesh entity to create an edge from.


    .. cpp:function:: double length() const
    
        Compute Euclidean length of edge
        
        *Returns*
            double
                Euclidean length of edge.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(2, 2);
                Edge edge(mesh, 0);
                info("%g", edge.length());
        
            output::
        
                0.5


    .. cpp:function:: double dot(const Edge& edge) const
    
        Compute dot product between edge and other edge
        
        *Arguments*
            edge (:cpp:class:`Edge`)
                Another edge.
        
        *Returns*
            double
                The dot product.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(2, 2);
                Edge edge1(mesh, 0);
                Edge edge2(mesh, 1);
                info("%g", edge1.dot(edge2));
        
            output::
        
                0.25


.. cpp:class:: EdgeFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    An EdgeFunction is a :cpp:class:`MeshFunction` of topological dimension 1.


