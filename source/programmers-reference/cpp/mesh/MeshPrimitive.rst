.. Documentation for the header file dolfin/mesh/MeshPrimitive.h

.. _programmers_reference_cpp_mesh_Mesh:

MeshPrimitive.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: //AABB_tree. Uses conversion operator in dolfin::Point to create a certain
                                                                                                 //CGAL Point_3 type.
                                                                                                 Point_3 reference_point() const
    
        Provides a reference point required by the Primitive concept of CGAL

    .. cpp:function:: //reasons. Explanation: We use a modified AABB_tree, in which the local BBox
                                                                                                   //functor class has been redefined to use the bbox function of dolfin mesh entities.
                                                                                                   //been used, which means that we would have had to convert dolfin cells into
                                                                                                   //CGAL primitives only to initialize the tree, which is probably very costly
                                                                                                   //for 1 million of triangles.
                                                                                                   //    CGAL::Bbox_3 bbox () const
    
        *Not* required by the CGAL primitive concept, but added for efficieny
        Otherwise the bbox function of the Datum object (see below) would have

