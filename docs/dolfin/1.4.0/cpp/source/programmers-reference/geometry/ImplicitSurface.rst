
.. Documentation for the header file dolfin/geometry/ImplicitSurface.h

.. _programmers_reference_cpp_geometry_implicitsurface:

ImplicitSurface.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ImplicitSurface

    This class is used to define a surface via a function f(x) -> R,
    where for a point y on the surface f(y) = 0.
    WARNING: This class is experimental and likely to change.


    .. cpp:function:: ImplicitSurface(Sphere s, std::string type)
    
        Create an isosurface
        
        *Arguments*
            s (Sphere)
                Bounding sphere for surface.
        
            type (std::string)
                Isosurface type. One of "manifold", "manifold_with_boundary"
                or "non_manifold".
        
        *Example*
            .. code-block:: c++
        
                IsoSurface isosurface(Sphere(Point(0.0, 0.2, 0.4), 1.0),
                                      "manifold_with_boundary");
        


    .. cpp:function:: double f0(const Point& point) const = 0
    
        Signed distance function surface. If f0(p) = 0, the point p is
        possibly on the surface, which case ImplicitSurface::f1 can be
        called to check.


    .. cpp:function:: double f1(const Point& point) const
    
        For a point for which f0 \approx 0, return <= is point is on
        is on the surface.  Can be used for creating open surfaces by
        discarding with any artificial closure.


