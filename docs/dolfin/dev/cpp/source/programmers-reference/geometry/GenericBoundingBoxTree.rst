
.. Documentation for the header file dolfin/geometry/GenericBoundingBoxTree.h

.. _programmers_reference_cpp_geometry_genericboundingboxtree:

GenericBoundingBoxTree.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericBoundingBoxTree

    Base class for bounding box implementations (envelope-letter
    design)


    .. cpp:function:: GenericBoundingBoxTree()
    
        Constructor


    .. cpp:function:: void build(const Mesh& mesh, std::size_t tdim)
    
        Build bounding box tree for mesh entites of given dimension


    .. cpp:function:: void build(const std::vector<Point>& points)
    
        Build bounding box tree for point cloud


    .. cpp:function:: std::vector<unsigned int> compute_collisions(const Point& point) const
    
        Compute all collisions between bounding boxes and given :cpp:class:`Point`


    .. cpp:function:: std::vector<unsigned int> compute_entity_collisions(const Point& point, const Mesh& mesh) const
    
        Compute all collisions between entities and given :cpp:class:`Point`


    .. cpp:function:: unsigned int compute_first_collision(const Point& point) const
    
        Compute first collision between bounding boxes and given :cpp:class:`Point`


    .. cpp:function:: unsigned int compute_first_entity_collision(const Point& point, const Mesh& mesh) const
    
        Compute first collision between entities and given :cpp:class:`Point`


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_entity(const Point& point, const Mesh& mesh) const
    
        Compute closest entity and distance to given :cpp:class:`Point`


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_point(const Point& point) const
    
        Compute closest point and distance to given :cpp:class:`Point`


    .. cpp:function:: void compute_collisions(const Point& point, unsigned int node, std::vector<unsigned int>& entities) const
    
        Compute collisions (recursive)


    .. cpp:function:: void compute_entity_collisions(const Point& point, unsigned int node, std::vector<unsigned int>& entities, const Mesh& mesh) const
    
        Compute entity collisions (recursive)


    .. cpp:function:: unsigned int compute_first_collision(const Point& point, unsigned int node) const
    
        Compute first collision (recursive)


    .. cpp:function:: unsigned int compute_first_entity_collision(const Point& point, unsigned int node, const Mesh& mesh) const
    
        Compute first entity collision (recursive)


    .. cpp:function:: void compute_closest_entity(const Point& point, unsigned int node, const Mesh& mesh, unsigned int& closest_entity, double& R2) const
    
        Compute closest entity (recursive)


    .. cpp:function:: void compute_closest_point(const Point& point, unsigned int node, unsigned int& closest_point, double& R2) const
    
        Compute closest point (recursive)


