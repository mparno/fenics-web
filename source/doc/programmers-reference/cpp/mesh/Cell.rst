.. Documentation for the header file dolfin/mesh/Cell.h

.. _programmers_reference_cpp_mesh_cell:

Cell.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Cell

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
    A Cell is a MeshEntity of topological codimension 0.

    .. cpp:function:: Cell()
    
        Create empty cell

    .. cpp:function:: Cell(const Mesh& mesh, uint index)
    
        Create cell on given mesh with given index

    .. cpp:function:: CellType::Type type() const
    
        Return type of cell

    .. cpp:function:: double orientation() const
    
        Compute orientation of cell (0 is right, 1 is left)

    .. cpp:function:: double volume() const
    
        Compute (generalized) volume of cell

    .. cpp:function:: double diameter() const
    
        Compute diameter of cell

    .. cpp:function:: double normal(uint facet, uint i) const
    
        Compute component i of normal of given facet with respect to the cell

    .. cpp:function:: Point normal(uint facet) const
    
        Compute normal of given facet with respect to the cell

    .. cpp:function:: double facet_area(uint facet) const
    
        Compute the area/length of given facet with respect to the cell

    .. cpp:function:: void order(const MeshFunction<uint>* global_vertex_indices)
    
        Order entities locally

    .. cpp:function:: bool ordered(const MeshFunction<uint>* global_vertex_indices) const
    
        Check if entities are ordered

.. cpp:class:: CellIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
    A CellIterator is a MeshEntityIterator of topological codimension 0.

.. cpp:class:: CellFunction

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
    A CellFunction is a MeshFunction of topological codimension 0.

