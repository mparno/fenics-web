
.. Documentation for the header file dolfin/generation/CSGPrimitives2D.h

.. _programmers_reference_cpp_generation_csgprimitives2d:

CSGPrimitives2D.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CSGPrimitive2D

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive`
        
    Base class for 2D primitives


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of geometry


.. cpp:class:: Circle

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive2D`
        
    This class describes a 2D circle which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Circle(double x0, double x1, double r, std::size_t fragments=32)
    
        Create circle at x = (x0, x1) with radius r.
        
        *Arguments*
            x0 (double)
                x0-coordinate of center.
            x1 (double)
                x1-coordinate of center.
            r (double)
                radius.
            fragments (std::size_t)
                number of fragments.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


    .. cpp:function:: Point center() const
    
        Return center of circle


    .. cpp:function:: double radius() const
    
        Return radius of circle


    .. cpp:function:: std::size_t fragments() const
    
        Return number of fragments around the circle


.. cpp:class:: Ellipse

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive2D`
        
    This class describes a 2D ellipse which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Ellipse(double x0, double x1, double a, double b, std::size_t fragments=32)
    
        Create ellipse at x = (x0, x1) with horizontal semi-axis a and
        vertical semi-axis b.
        
        *Arguments*
            x0 (double)
                x0-coordinate of center.
            x1 (double)
                x1-coordinate of center.
            a (double)
                horizontal semi-axis.
            b (double)
                vertical semi-axis.
            fragments (std::size_t)
                number of fragments.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


    .. cpp:function:: Point center() const
    
        Return center of ellipse


    .. cpp:function:: double a() const
    
        Return horizontal semi-axis


    .. cpp:function:: double b() const
    
        Return vertical semi-axis


    .. cpp:function:: std::size_t fragments() const
    
        Return number of fragments around the ellipse


.. cpp:class:: Rectangle

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive2D`
        
    This class describes a 2D rectangle which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Rectangle(double x0, double x1, double y0, double y1)
    
        Create rectangle defined by two opposite corners
        x = (x0, x1) and y = (y0, y1).
        
        *Arguments*
            x0 (double)
                x0-coordinate of first corner.
            x1 (double)
                x1-coordinate of first corner.
            y0 (double)
                y0-coordinate of second corner.
            y1 (double)
                y1-coordinate of second corner.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


    .. cpp:function:: Point first_corner() const
    
        Return first corner


    .. cpp:function:: Point second_corner() const
    
        Return second corner


.. cpp:class:: Polygon

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive2D`
        
    This class describes a 2D polygon which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


    .. cpp:function:: const std::vector<Point>& vertices() const
    
        Return vertices in polygon


