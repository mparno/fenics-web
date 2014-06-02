
.. Documentation for the header file dolfin/mesh/CellType.h

.. _programmers_reference_cpp_mesh_celltype:

CellType.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CellType

    This class provides a common interface for different cell types.
    Each cell type implements mesh functionality that is specific to
    a certain type of cell.


    .. cpp:function:: CellType(Type cell_type, Type facet_type)
    
        Constructor


    .. cpp:function:: static CellType* create(Type type)
    
        Create cell type from type (factory function)


    .. cpp:function:: static CellType* create(std::string type)
    
        Create cell type from string (factory function)


    .. cpp:function:: static Type string2type(std::string type)
    
        Convert from string to cell type


    .. cpp:function:: static std::string type2string(Type type)
    
        Convert from cell type to string


    .. cpp:function:: Type cell_type() const
    
        Return type of cell


    .. cpp:function:: Type facet_type() const
    
        Return type of cell for facets


    .. cpp:function:: std::size_t dim() const = 0
    
        Return topological dimension of cell


    .. cpp:function:: std::size_t num_entities(std::size_t dim) const = 0
    
        Return number of entitites of given topological dimension


    .. cpp:function:: std::size_t num_vertices(std::size_t dim) const = 0
    
        Return number of vertices for entity of given topological dimension


    .. cpp:function:: std::size_t orientation(const Cell& cell) const = 0
    
        Return orientation of the cell (assuming flat space)


    .. cpp:function:: std::size_t orientation(const Cell& cell, const Point& up) const
    
        Return orientation of the cell relative to given up direction


    .. cpp:function:: void create_entities(std::vector<std::vector<unsigned int> >& e, std::size_t dim, const unsigned int* v) const = 0
    
        Create entities e of given topological dimension from vertices v


    .. cpp:function:: void refine_cell(Cell& cell, MeshEditor& editor, std::size_t& current_cell) const = 0
    
        Refine cell uniformly


    .. cpp:function:: double volume(const MeshEntity& entity) const = 0
    
        Compute (generalized) volume of mesh entity


    .. cpp:function:: double diameter(const MeshEntity& entity) const = 0
    
        Compute diameter of mesh entity


    .. cpp:function:: double inradius(const Cell& cell) const
    
        Compute inradius of cell


    .. cpp:function:: double radius_ratio(const Cell& cell) const
    
        Compute dim*inradius/circumradius for given cell


    .. cpp:function:: double squared_distance(const Cell& cell, const Point& point) const = 0
    
        Compute squared distance to given point


    .. cpp:function:: double normal(const Cell& cell, std::size_t facet, std::size_t i) const = 0
    
        Compute component i of normal of given facet with respect to the cell


    .. cpp:function:: Point normal(const Cell& cell, std::size_t facet) const = 0
    
        Compute of given facet with respect to the cell


    .. cpp:function:: Point cell_normal(const Cell& cell) const = 0
    
        Compute normal to given cell (viewed as embedded in 3D)


    .. cpp:function:: double facet_area(const Cell& cell, std::size_t facet) const = 0
    
        Compute the area/length of given facet with respect to the cell


    .. cpp:function:: void order(Cell& cell, const std::vector<std::size_t>& local_to_global_vertex_indices) const = 0
    
        Order entities locally


    .. cpp:function:: bool ordered(const Cell& cell, const std::vector<std::size_t>& local_to_global_vertex_indices) const
    
        Check if entities are ordered


    .. cpp:function:: bool collides(const Cell& cell, const Point& point) const = 0
    
        Check whether given point collides with cell


    .. cpp:function:: bool collides(const Cell& cell, const MeshEntity& entity) const = 0
    
        Check whether given entity collides with cell


    .. cpp:function:: std::vector<double> triangulate_intersection(const Cell& c0, const Cell& c1) const = 0
    
        Compute triangulation of intersection of two cells


    .. cpp:function:: std::string description(bool plural) const = 0
    
        Return description of cell type


