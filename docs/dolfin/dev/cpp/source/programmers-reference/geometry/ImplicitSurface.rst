
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
        


    .. cpp:function:: double operator()(const Point& point) const = 0
    
        Return value of isosurfacce function. This function is
        overloaed by the user. For closed volumes, the function should
        evaluate to < 0 on the interior.
        
        *Arguments*
            point (Point)
                The point at which to evaluate the isosurface function.
                or "non_manifold".
        
        *Returns*
            double
                Isosurface function value.


    .. cpp:function:: bool on_surface(const Point& point) const
    
        Test whether a point for which f \approx 0 is on the
        surface. Used for creating open surfaces by discarding with
        any artificial closure.


