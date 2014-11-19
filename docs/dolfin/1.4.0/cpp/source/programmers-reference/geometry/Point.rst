
.. Documentation for the header file dolfin/geometry/Point.h

.. _programmers_reference_cpp_geometry_point:

Point.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Point

    A Point represents a point in :math:`\mathbb{R}^3` with
    coordinates :math:`x, y, z,` or alternatively, a vector in
    :math:`\mathbb{R}^3`, supporting standard operations like the
    norm, distances, scalar and vector products etc.


    .. cpp:function:: Point(const double x=0.0, const double y=0.0, const double z=0.0)
    
        Create a point at (x, y, z). Default value (0, 0, 0).
        
        *Arguments*
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.
            z (double)
                The z-coordinate.


    .. cpp:function:: Point(std::size_t dim, const double* x)
    
        Create point from array
        
        *Arguments*
            dim (std::size_t)
                Dimension of the array.
            x (double)
                The array to create a Point from.


    .. cpp:function:: Point(const Point& p)
    
        Copy constructor
        
        *Arguments*
            p (:cpp:class:`Point`)
                The object to be copied.


    .. cpp:function:: double& operator[] (std::size_t i)
    
        Return address of coordinate in direction i
        
        *Arguments*
            i (std::size_t)
                Direction.
        
        *Returns*
            double
                Adress of coordinate in the given direction.


    .. cpp:function:: double operator[] (std::size_t i) const
    
        Return coordinate in direction i
        
        *Arguments*
            i (std::size_t)
                Direction.
        
        *Returns*
            double
                The coordinate in the given direction.


    .. cpp:function:: double x() const
    
        Return x-coordinate
        
        *Returns*
            double
                The x-coordinate.


    .. cpp:function:: double y() const
    
        Return y-coordinate
        
        *Returns*
            double
                The y-coordinate.


    .. cpp:function:: double z() const
    
        Return z-coordinate
        
        *Returns*
            double
                The z-coordinate.


    .. cpp:function:: double* coordinates()
    
        Return coordinate array
        
        *Returns*
            list of doubles
                The coordinates.


    .. cpp:function:: const double* coordinates() const
    
        Return coordinate array (const. version)
        
        *Returns*
            list of doubles
                The coordinates.


    .. cpp:function:: Point operator+ (const Point& p) const
    
        Compute sum of two points


    .. cpp:function:: Point operator- (const Point& p) const
    
        Compute difference of two points


    .. cpp:function:: const Point& operator+= (const Point& p)
    
        Add given point


    .. cpp:function:: const Point& operator-= (const Point& p)
    
        Subtract given point


    .. cpp:function:: Point operator- ()
    
        Unary minus


    .. cpp:function:: Point operator* (double a) const
    
        Multiplication with scalar


    .. cpp:function:: const Point& operator*= (double a)
    
        Incremental multiplication with scalar


    .. cpp:function:: Point operator/ (double a) const
    
        Division by scalar


    .. cpp:function:: const Point& operator/= (double a)
    
        Incremental division by scalar


    .. cpp:function:: const Point& operator= (const Point& p)
    
        Assignment operator


    .. cpp:function:: double squared_distance(const Point& p) const
    
        Compute squared distance to given point
        
        *Arguments*
            p (:cpp:class:`Point`)
                The point to compute distance to.
        
        *Returns*
            double
                The squared distance.
        
        *Example*
            .. code-block:: c++
        
                Point p1(0, 4, 0);
                Point p2(2, 0, 4);
                info("%g", p1.squared_distance(p2));
        
            output::
        
                6


    .. cpp:function:: double distance(const Point& p) const
    
        Compute distance to given point
        
        *Arguments*
            p (:cpp:class:`Point`)
                The point to compute distance to.
        
        *Returns*
            double
                The distance.
        
        *Example*
            .. code-block:: c++
        
                Point p1(0, 4, 0);
                Point p2(2, 0, 4);
                info("%g", p1.distance(p2));
        
            output::
        
                6


    .. cpp:function:: double norm() const
    
        Compute norm of point representing a vector from the origin
        
        *Returns*
            double
                The (Euclidean) norm of the vector from the origin to
                the point.
        
        *Example*
            .. code-block:: c++
        
                Point p(1.0, 2.0, 2.0);
                info("%g", p.norm());
        
            output::
        
                3


    .. cpp:function:: double squared_norm() const
    
        Compute norm of point representing a vector from the origin
        
        *Returns*
            double
                The squared (Euclidean) norm of the vector from the
                origin of the point.
        
        *Example*
            .. code-block:: c++
        
                Point p(1.0, 2.0, 2.0);
                info("%g", p.squared_norm());
        
            output::
        
                9


    .. cpp:function:: const Point cross(const Point& p) const
    
        Compute cross product with given vector
        
        *Arguments*
            p (:cpp:class:`Point`)
                Another point.
        
        *Returns*
            Point
                The cross product.


    .. cpp:function:: double dot(const Point& p) const
    
        Compute dot product with given vector
        
        *Arguments*
            p (:cpp:class:`Point`)
                Another point.
        
        *Returns*
            double
                The dot product.
        
        *Example*
            .. code-block:: c++
        
                Point p1(1.0, 4.0, 8.0);
                Point p2(2.0, 0.0, 0.0);
                info("%g", p1.dot(p2));
        
            output::
        
                2


    .. cpp:function:: Point rotate(const Point& a, double theta) const
    
        Rotate around a given axis
        
        *Arguments*
            a (:cpp:class:`Point`)
                The axis to rotate around. Must be unit length.
            theta (_double_)
                The rotation angle.
        
        *Returns*
            Point
                The rotated point.


    .. cpp:function:: std::string str(bool verbose=false) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


    .. cpp:function:: Point operator*(double a, const Point& p)
    
        Multiplication with scalar


