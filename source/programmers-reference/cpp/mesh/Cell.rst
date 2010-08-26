.. Documentation for the header file dolfin/mesh/Cell.h

.. _programmers_reference_cpp_mesh_cell:

Cell.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Cell

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
        A Cell is a MeshEntity of topological codimension 0.

    .. cpp:function:: Cell() : MeshEntity()
    
        Create empty cell

    .. cpp:function:: Cell(const Mesh& mesh, uint index) : MeshEntity(mesh, mesh.topology().dim(), index)
    
        Create cell on given mesh with given index

    .. cpp:function:: inline CellType::Type type() const
    
        Return type of cell

    .. cpp:function:: inline Point normal(uint facet) const
    
        Compute normal of given facet with respect to the cell

    .. cpp:function:: inline bool ordered(MeshFunction<uint>* global_vertex_indices)
    
        Check if entities are ordered

    .. cpp:function:: inline double diameter() const
    
        Compute diameter of cell

    .. cpp:function:: inline double facet_area(uint facet) const
    
        Compute the area/length of given facet with respect to the cell

    .. cpp:function:: inline double normal(uint facet, uint i) const
    
        Compute component i of normal of given facet with respect to the cell

    .. cpp:function:: inline double orientation() const
    
        Compute orientation of cell (0 is right, 1 is left)

    .. cpp:function:: inline double volume() const
    
        Compute (generalized) volume of cell

    .. cpp:function:: inline void order(MeshFunction<uint>* global_vertex_indices)
    
        Order entities locally

    .. cpp:function:: ~Cell()
    
        Destructor

.. cpp:class:: CellIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
        A CellIterator is a MeshEntityIterator of topological codimension 0.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
        A CellFunction is a MeshFunction of topological codimension 0.

