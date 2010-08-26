.. Documentation for the header file dolfin/mesh/PrimitiveTraits.h

.. _programmers_reference_cpp_mesh_primitivetraits:

PrimitiveTraits.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: template <typename Primitive_, typename Kernel> struct PrimitiveTraits
    
        Forward declaration for a general Traits class. This traits class is
        supposed to provide a datum function, which returns a geometric primitive
        object, which type corresponds to the primitive type (Point, PointCell,
        Tetrahedron(Cell) etc.) and the passed geometric CGAL kernel.

