
.. Documentation for the header file dolfin/generation/CSGPrimitives3D.h

.. _programmers_reference_cpp_generation_csgprimitives3d:

CSGPrimitives3D.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CSGPrimitive3D

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive`
        
    Base class for 3D primitives


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of geometry


.. cpp:class:: Sphere

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive3D`
        
    This class describes a 3D sphere which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Sphere(Point center, double radius, std::size_t slices=16)
    
        Create sphere at x = (x0, x1, x2) with radius r.
        
        *Arguments*
            x0 (double)
                x0-coordinate of center.
            x1 (double)
                x1-coordinate of center.
            x2 (double)
                x2-coordinate of center.
            r (double)
                radius.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: Box

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive3D`
        
    This class describes a 3D box which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Box(double x0, double x1, double x2, double y0, double y1, double y2)
    
        Create box defined by two opposite corners
        x = (x0, x1, x2) and y = (y0, y1, y2).
        
        *Arguments*
            x0 (double)
                x0-coordinate of first corner.
            x1 (double)
                x1-coordinate of first corner.
            x2 (double)
                x2-coordinate of first corner.
            y0 (double)
                y0-coordinate of second corner.
            y1 (double)
                y1-coordinate of second corner.
            y2 (double)
                y2-coordinate of second corner.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: Cone

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive3D`
        
    This class describes a 3D cone which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Cone(Point top, Point bottom, double top_radius, double bottom_radius, std::size_t slices=32)
    
        Create cone defined by upper and lower center
        and radius respectively.
        
        *Arguments*
            top (Point)
                Center at top of cone.
            top_radius(double)
                Radius bottom of cone.
            bottom(Point)
                Center at top of cone.
            bottom_radius (double)
                radius at top of cone.
            slices (std::size_t)
                number of faces on the side when generating a
                polyhedral approximation.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: Cylinder

    *Parent class(es)*
    
        * :cpp:class:`Cone`
        
    This class describes a 3D cylinder which can be used to build
    geometries using Constructive Solid Geometry (CSG). A cylinder
    is here just a special case of a cone.


    .. cpp:function:: Cylinder(Point top, Point bottom, double r, std::size_t slices=32)
    
        Create cylinder defined by upper and lower center
        and radius respectively.
        
        *Arguments*
            top (Point)
                Center at top of cylinder.
            bottom(Point)
                Center at top of cylinder.
            r (double)
                radius of cylinder.
            slices (std::size_t)
                number of faces on the side when generating a
                polyhedral approximation.


.. cpp:class:: Tetrahedron

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive3D`
        
    This class describes a Tetrahedron which can be used to build
    geometries using Constructive Solid Geometry (CSG).


    .. cpp:function:: Tetrahedron(Point x0, Point x1, Point x2, Point x3)
    
        Create tetrahedron defined by four corner points.
        
        *Arguments*
            x0 (Point)
                Point.
            x1 (Point)
                Point.
            x2 (Point)
                Point.
            x3 (Point)
                Point.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


.. cpp:class:: Surface3D

    *Parent class(es)*
    
        * :cpp:class:`CSGPrimitive3D`
        
    This class describes a 3D surface loaded from file.
    The supported file types


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation


