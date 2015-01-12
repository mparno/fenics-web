
.. Documentation for the header file dolfin/geometry/MeshPointIntersection.h

.. _programmers_reference_cpp_geometry_meshpointintersection:

MeshPointIntersection.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshPointIntersection

    This class represents an intersection between a :cpp:class:`Mesh` and a
    :cpp:class:`Point`. The resulting intersection is stored as a list of zero
    or more cells.


    .. cpp:function:: MeshPointIntersection(const Mesh& mesh, const Point& point)
    
        Compute intersection between mesh and point


    .. cpp:function:: const std::vector<unsigned int>& intersected_cells() const
    
        Return the list of (local) indices for intersected cells


