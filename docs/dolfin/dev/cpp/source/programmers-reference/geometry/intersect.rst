
.. Documentation for the header file dolfin/geometry/intersect.h

.. _programmers_reference_cpp_geometry_intersect:

intersect.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: boost::shared_ptr<const MeshPointIntersection> intersect(const Mesh& mesh, const Point& point)

    Compute and return intersection between :cpp:class:`Mesh` and :cpp:class:`Point`.
    
    *Arguments*
        mesh (:cpp:class:`Mesh`)
            The mesh to be intersected.
        point (:cpp:class:`Point`)
            The point to be intersected.
    
    *Returns*
        :cpp:class:`MeshPointIntersection`
            The intersection data.


