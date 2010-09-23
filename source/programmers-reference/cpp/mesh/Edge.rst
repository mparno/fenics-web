.. Documentation for the header file dolfin/mesh/Edge.h

.. _programmers_reference_cpp_mesh_edge:

Edge.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Edge

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
    An Edge is a MeshEntity of topological dimension 1.

    .. cpp:function:: Edge(const Mesh& mesh, uint index)
    
        Create edge on given mesh

    .. cpp:function:: Edge(MeshEntity& entity)
    
        Create edge from mesh entity

    .. cpp:function:: double length()
    
        Compute Euclidean length of edge

.. cpp:class:: EdgeIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
    An EdgeIterator is a MeshEntityIterator of topological dimension 1.

.. cpp:class:: EdgeFunction

    *Parent class*
    
        * :cpp:class:`MeshFunction`
        
    An EdgeFunction is a MeshFunction of topological dimension 1.

