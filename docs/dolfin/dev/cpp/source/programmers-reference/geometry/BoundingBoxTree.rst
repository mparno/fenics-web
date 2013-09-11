
.. Documentation for the header file dolfin/geometry/BoundingBoxTree.h

.. _programmers_reference_cpp_geometry_boundingboxtree:

BoundingBoxTree.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BoundingBoxTree

    This class implements a (distributed) axis aligned bounding box
    tree (AABB tree). Bounding box trees can be created from meshes
    and [other data structures, to be filled in].


    .. cpp:function:: BoundingBoxTree()
    
        Create empty bounding box tree


    .. cpp:function:: void build(const Mesh& mesh)
    
        Build bounding box tree for cells of mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh for which to compute the bounding box tree.


    .. cpp:function:: void build(const Mesh& mesh, std::size_t tdim)
    
        Build bounding box tree for mesh entities of given dimension.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh for which to compute the bounding box tree.
            dimension (std::size_t)
                The entity dimension (topological dimension) for which
                to compute the bounding box tree.


    .. cpp:function:: void build(const std::vector<Point>& points, std::size_t gdim)
    
        Build bounding box tree for point cloud.
        
        *Arguments*
            points (std::vector<:cpp:class:`Point`>)
                The list of points.
            gdim (std::size_t)
                The geometric dimension.


    .. cpp:function:: std::vector<unsigned int> compute_collisions(const Point& point) const
    
        Compute all collisions between bounding boxes and given :cpp:class:`Point`.
        
        *Returns*
            std::vector<unsigned int>
                A list of local indices for entities contained in
                (leaf) bounding boxes that collide with (intersect)
                the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.


    .. cpp:function:: std::vector<unsigned int> compute_entity_collisions(const Point& point, const Mesh& mesh) const
    
        Compute all collisions between entities and given :cpp:class:`Point`.
        
        *Returns*
            std::vector<unsigned int>
                A list of local indices for entities that collide with
                (intersect) the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: unsigned int compute_first_collision(const Point& point) const
    
        Compute first collision between bounding boxes and given :cpp:class:`Point`.
        
        *Returns*
            unsigned int
                The local index for the first found entity contained
                in a (leaf) bounding box that collides with
                (intersects) the given point. If not found,
                std::numeric_limits<unsigned int>::max() is returned.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.


    .. cpp:function:: unsigned int compute_first_entity_collision(const Point& point, const Mesh& mesh) const
    
        Compute first collision between entities and given :cpp:class:`Point`.
        
        *Returns*
            unsigned int
                The local index for the first found entity that
                collides with (intersects) the given point. If not
                found, std::numeric_limits<unsigned int>::max() is
                returned.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_entity(const Point& point, const Mesh& mesh) const
    
        Compute closest entity to given :cpp:class:`Point`.
        
        *Returns*
            unsigned int
                The local index for the entity that is closest to the
                point. If more than one entity is at the same distance
                (or point contained in entity), then the first entity
                is returned.
            double
                The distance to the closest entity.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_point(const Point& point) const
    
        Compute closest point to given :cpp:class:`Point`. This function assumes
        that the tree has been built for a point cloud.
        
        Developer note: This function should not be confused with
        computing the closest point in all entities of a mesh. That
        function could be added with relative ease since we actually
        compute the closest points to get the distance in the above
        function (compute_closest_entity) inside the specialized
        implementations in TetrahedronCell.cpp etc.
        
        *Returns*
            unsigned int
                The local index for the point that is closest to the
                point. If more than one point is at the same distance
                (or point contained in entity), then the first point
                is returned.
            double
                The distance to the closest point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.


