.. Documentation for the header file dolfin/mesh/PrimitiveIntersector.h

.. _programmers_reference_cpp_mesh_primitiveintersector:

PrimitiveIntersector.h
======================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: PrimitiveIntersector

    This class implements an intersection detection, detecting
    whether two given (arbitrary) meshentities intersect.

    .. cpp:function:: static bool do_intersect(const MeshEntity& entity_1, const MeshEntity& entity_2)
    
        Computes whether two mesh entities intersect using an inexact
        geometry kernel which is fast but may suffer from floating
        point precision

    .. cpp:function:: static bool do_intersect(const MeshEntity& entity_1, const Point& point)
    
        Computes whether a mesh entity and point intersect using an
        inexact geometry kernel which is fast but may suffer from
        floating point precision

    .. cpp:function:: static bool do_intersect_exact(const MeshEntity& entity_1, const MeshEntity& entity_2)
    
        Computes whether two mesh entities intersect using an exact
        geometry kernel which is slow but always correct

    .. cpp:function:: static bool do_intersect_exact(const MeshEntity& entity_1, const Point& point)
    
        Computes whether a mesh entity and point intersect using an
        exact geometry kernel which is slow but always correct

