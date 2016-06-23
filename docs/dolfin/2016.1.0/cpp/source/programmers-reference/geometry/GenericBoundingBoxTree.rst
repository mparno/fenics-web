
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


    .. cpp:function:: static std::shared_ptr<GenericBoundingBoxTree> create(unsigned int dim)
    
        Factory function returning (empty) tree of appropriate dimension


    .. cpp:function:: void build(const Mesh& mesh, std::size_t tdim)
    
        Build bounding box tree for mesh entities of given dimension


    .. cpp:function:: void build(const std::vector<Point>& points)
    
        Build bounding box tree for point cloud


    .. cpp:function:: std::vector<unsigned int> compute_collisions(const Point& point) const
    
        Compute all collisions between bounding boxes and :cpp:class:`Point`


    .. cpp:function:: std::pair<std::vector<unsigned int>, std::vector<unsigned int> > compute_collisions(const GenericBoundingBoxTree& tree) const
    
        Compute all collisions between bounding boxes and :cpp:class:`BoundingBoxTree`


    .. cpp:function:: std::vector<unsigned int> compute_entity_collisions(const Point& point, const Mesh& mesh) const
    
        Compute all collisions between entities and :cpp:class:`Point`


    .. cpp:function:: std::vector<unsigned int> compute_process_collisions(const Point& point) const
    
        Compute all collisions between processes and :cpp:class:`Point`


    .. cpp:function:: std::pair<std::vector<unsigned int>, std::vector<unsigned int> > compute_entity_collisions(const GenericBoundingBoxTree& tree, const Mesh& mesh_A, const Mesh& mesh_B) const
    
        Compute all collisions between entities and :cpp:class:`BoundingBoxTree`


    .. cpp:function:: unsigned int compute_first_collision(const Point& point) const
    
        Compute first collision between bounding boxes and :cpp:class:`Point`


    .. cpp:function:: unsigned int compute_first_entity_collision(const Point& point, const Mesh& mesh) const
    
        Compute first collision between entities and :cpp:class:`Point`


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_entity(const Point& point, const Mesh& mesh) const
    
        Compute closest entity and distance to :cpp:class:`Point`


    .. cpp:function:: std::pair<unsigned int, double> compute_closest_point(const Point& point) const
    
        Compute closest point and distance to :cpp:class:`Point`


    .. cpp:function:: std::string str(bool verbose=false)
    
        Print out for debugging


    .. cpp:function:: static void _compute_collisions(const GenericBoundingBoxTree& tree, const Point& point, unsigned int node, std::vector<unsigned int>& entities, const Mesh* mesh)
    
        Compute collisions with point (recursive)


    .. cpp:function:: static void _compute_collisions(const GenericBoundingBoxTree& A, const GenericBoundingBoxTree& B, unsigned int node_A, unsigned int node_B, std::vector<unsigned int>& entities_A, std::vector<unsigned int>& entities_B, const Mesh* mesh_A, const Mesh* mesh_B)
    
        Compute collisions with tree (recursive)


    .. cpp:function:: static unsigned int _compute_first_collision(const GenericBoundingBoxTree& tree, const Point& point, unsigned int node)
    
        Compute first collision (recursive)


    .. cpp:function:: static unsigned int _compute_first_entity_collision(const GenericBoundingBoxTree& tree, const Point& point, unsigned int node, const Mesh& mesh)
    
        Compute first entity collision (recursive)


    .. cpp:function:: static void _compute_closest_entity(const GenericBoundingBoxTree& tree, const Point& point, unsigned int node, const Mesh& mesh, unsigned int& closest_entity, double& R2)
    
        Compute closest entity (recursive)


    .. cpp:function:: static void _compute_closest_point(const GenericBoundingBoxTree& tree, const Point& point, unsigned int node, unsigned int& closest_point, double& R2)
    
        Compute closest point (recursive)


