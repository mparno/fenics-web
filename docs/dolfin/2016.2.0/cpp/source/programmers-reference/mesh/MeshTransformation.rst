
.. Documentation for the header file dolfin/mesh/MeshTransformation.h

.. _programmers_reference_cpp_mesh_meshtransformation:

MeshTransformation.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: static void translate(Mesh& mesh, const Point& point)

    Translate mesh according to a given vector.
    
    *Arguments*
        mesh (:cpp:class:`Mesh`)
            The mesh
        point (Point)
            The vector defining the translation.


.. cpp:function:: static void rescale(Mesh& mesh, const double scale, const Point& center)

    Rescale mesh by a given scaling factor with respect to a center
    point.
    
    *Arguments*
        mesh (:cpp:class:`Mesh`)
            The mesh
         scale (double)
            The scaling factor.
         center (Point)
            The center of the scaling.


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


