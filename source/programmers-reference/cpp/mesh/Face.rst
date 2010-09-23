.. Documentation for the header file dolfin/mesh/Face.h

.. _programmers_reference_cpp_mesh_face:

Face.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Face

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
    A Face is a MeshEntity of topological dimension 2.

    .. cpp:function:: Face(const Mesh& mesh, uint index)
    
        Constructor

    .. cpp:function:: double area() const
    
        Calculate the area of the face (triangle)

.. cpp:class:: FaceIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
    A FaceIterator is a MeshEntityIterator of topological dimension 2.

.. cpp:class:: FaceFunction

    *Parent class*
    
        * :cpp:class:`MeshFunction`
        
    A FaceFunction is a MeshFunction of topological dimension 2.

