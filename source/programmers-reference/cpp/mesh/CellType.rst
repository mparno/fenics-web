.. Documentation for the header file dolfin/mesh/CellType.h

.. _programmers_reference_cpp_mesh_celltype:

CellType.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: CellType

    This class provides a common interface for different cell types.
    Each cell type implements mesh functionality that is specific to
    a certain type of cell.

    .. cpp:function:: enum Type
    
        Enum for different cell types

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

    .. cpp:function:: uint dim() const = 0
    
        Return topological dimension of cell

    .. cpp:function:: uint num_entities(uint dim) const = 0
    
        Return number of entitites of given topological dimension

    .. cpp:function:: uint num_vertices(uint dim) const = 0
    
        Return number of vertices for entity of given topological dimension

    .. cpp:function:: uint orientation(const Cell& cell) const = 0
    
        Return orientation of the cell

    .. cpp:function:: void create_entities(uint** e, uint dim, const uint* v) const = 0
    
        Create entities e of given topological dimension from vertices v

    .. cpp:function:: void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const = 0
    
        Refine cell uniformly

    .. cpp:function:: double volume(const MeshEntity& entity) const = 0
    
        Compute (generalized) volume of mesh entity

    .. cpp:function:: double diameter(const MeshEntity& entity) const = 0
    
        Compute diameter of mesh entity

    .. cpp:function:: double normal(const Cell& cell, uint facet, uint i) const = 0
    
        Compute component i of normal of given facet with respect to the cell

    .. cpp:function:: Point normal(const Cell& cell, uint facet) const = 0
    
        Compute of given facet with respect to the cell

    .. cpp:function:: double facet_area(const Cell& cell, uint facet) const = 0
    
        Compute the area/length of given facet with respect to the cell

    .. cpp:function:: void order(Cell& cell, const MeshFunction<uint>* global_vertex_indices) const = 0
    
        Order entities locally

    .. cpp:function:: bool ordered(const Cell& cell, MeshFunction<uint>* global_vertex_indices) const
    
        Check if entities are ordered

    .. cpp:function:: std::string description(bool plural) const = 0
    
        Return description of cell type

