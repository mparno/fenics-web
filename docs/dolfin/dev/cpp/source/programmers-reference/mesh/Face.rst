
.. Documentation for the header file dolfin/mesh/Face.h

.. _programmers_reference_cpp_mesh_face:

Face.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Face

    *Parent class(es)*
    
        * :cpp:class:`MeshEntity`
        
    A Face is a MeshEntity of topological dimension 2.


    .. cpp:function:: Face(const Mesh& mesh, uint index)
    
        Constructor


    .. cpp:function:: double area() const
    
        Calculate the area of the face (triangle)


    .. cpp:function:: double normal(uint i) const
    
        Compute component i of the normal to the face


    .. cpp:function:: Point normal() const
    
        Compute normal to the face


.. cpp:class:: FaceIterator

    *Parent class(es)*
    
        * :cpp:class:`MeshEntityIterator`
        
    A FaceIterator is a MeshEntityIterator of topological dimension 2.


.. cpp:class:: FaceFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A FaceFunction is a MeshFunction of topological dimension 2.


