
.. Documentation for the header file dolfin/mesh/MeshTransformation.h

.. _programmers_reference_cpp_mesh_meshtransformation:

MeshTransformation.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshTransformation

    This class implements various transformations of the coordinates
    of a mesh.


    .. cpp:function:: static void translate(Mesh& mesh, const Point& point)
    
        Translate mesh according to a given vector.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh
            point (Point)
                The vector defining the translation.


    .. cpp:function:: static void rotate(Mesh& mesh, double angle, std::size_t axis)
    
        Rotate mesh around a coordinate axis through center of mass
        of all mesh vertices
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            angle (double)
                The number of degrees (0-360) of rotation.
            axis (std::size_t)
                The coordinate axis around which to rotate the mesh.


    .. cpp:function:: static void rotate(Mesh& mesh, double angle, std::size_t axis, const Point& p)
    
        Rotate mesh around a coordinate axis through a given point
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            angle (double)
                The number of degrees (0-360) of rotation.
            axis (std::size_t)
                The coordinate axis around which to rotate the mesh.
            point (:cpp:class:`Point`)
                The point around which to rotate the mesh.


