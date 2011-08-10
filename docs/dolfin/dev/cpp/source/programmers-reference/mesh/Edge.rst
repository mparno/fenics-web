
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
        
    An Edge is a MeshEntity of topological dimension 1.


    .. cpp:function:: Edge(const Mesh& mesh, uint index)
    
        Create edge on given mesh


    .. cpp:function:: Edge(MeshEntity& entity)
    
        Create edge from mesh entity


    .. cpp:function:: double length() const
    
        Compute Euclidean length of edge


    .. cpp:function:: double dot(const Edge& edge) const
    
        Compute dot product between edge and other edge


.. cpp:class:: EdgeIterator

    *Parent class(es)*
    
        * :cpp:class:`MeshEntityIterator`
        
    An EdgeIterator is a MeshEntityIterator of topological dimension 1.


.. cpp:class:: EdgeFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    An EdgeFunction is a MeshFunction of topological dimension 1.


