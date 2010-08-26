.. Documentation for the header file dolfin/mesh/Point.h

.. _programmers_reference_cpp_mesh_Mesh:

Point.h
=======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Point

        A Point represents a point in R^3 with coordinates x, y, z, or,
        alternatively, a vector in R^3, supporting standard operations
        like the norm, distances, scalar and vector products etc.

    .. cpp:function:: Point operator* (double a) const
    
        Multiplication with scalar

    .. cpp:function:: Point operator+ (const Point& p) const
    
        Compute sum of two points

    .. cpp:function:: Point operator- (const Point& p) const
    
        Compute difference of two points

    .. cpp:function:: Point operator/ (double a) const
    
        Division by scalar

    .. cpp:function:: Point(const Point& p)
    
        Copy constructor

    .. cpp:function:: Point(const double x = 0.0, const double y = 0.0, const double z =0.0)
    
        Create a point at (x, y, z)

    .. cpp:function:: Point(uint dim, const double* x)
    
        Create point from array

    .. cpp:function:: const Point cross(const Point& p) const
    
        Compute cross product with given vector

    .. cpp:function:: const Point& operator*= (double a)
    
        Incremental multiplication with scalar

    .. cpp:function:: const Point& operator+= (const Point& p)
    
        Add given point

    .. cpp:function:: const Point& operator-= (const Point& p)
    
        Subtract given point

    .. cpp:function:: const Point& operator/= (double a)
    
        Incremental division by scalar

    .. cpp:function:: const Point& operator= (const Point& p)
    
        Assignment operator

    .. cpp:function:: const double* coordinates() const
    
        Return coordinate array

    .. cpp:function:: double distance(const Point& p) const
    
        Compute distance to given point

    .. cpp:function:: double dot(const Point& p) const
    
        Compute dot product with given vector

    .. cpp:function:: double norm() const
    
        Compute norm of point representing a vector from the origin

    .. cpp:function:: double operator[] (uint i) const
    
        Return coordinate in direction i

    .. cpp:function:: double x() const
    
        Return x-coordinate

    .. cpp:function:: double y() const
    
        Return y-coordinate

    .. cpp:function:: double z() const
    
        Return z-coordinate

    .. cpp:function:: double& operator[] (uint i)
    
        Return address of coordinate in direction i

    .. cpp:function:: inline Point operator*(double a, const Point& p)
    
        Multiplication with scalar

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: template <typename Kernel>
                                                 CGAL::Bbox_3  bbox()
    
        Provides a CGAL bounding box, using conversion operator.

    .. cpp:function:: template <typename Kernel>
                                                 Point (const CGAL::Point_3<Kernel> & point)
    
        Constructor taking a CGAL::Point_3. Allows conversion from CGAL Point_3 class to Point class.

    .. cpp:function:: template <typename Kernel>
                                                 operator CGAL::Point_3<Kernel>() const
    
        Conversion operator to appropriate CGAL Point_3 class.

    .. cpp:function:: ~Point()
    
        Destructor

