.. Documentation for the header file dolfin/mesh/Edge.h

.. _programmers_reference_cpp_mesh_Mesh:

Edge.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Edge

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
        An Edge is a MeshEntity of topological dimension 1.

    .. cpp:function:: Edge(MeshEntity& entity) : MeshEntity(entity.mesh(), 1, entity.index())
    
        Create edge from mesh entity

    .. cpp:function:: Edge(const Mesh& mesh, uint index) : MeshEntity(mesh, 1, index)
    
        Create edge on given mesh

    .. cpp:function:: double length()
    
        Compute Euclidean length of edge

    .. cpp:function:: ~Edge()
    
        Destructor

.. cpp:class:: EdgeIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
        An EdgeIterator is a MeshEntityIterator of topological dimension 1.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
        An EdgeFunction is a MeshFunction of topological dimension 1.

