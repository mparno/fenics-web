
.. Documentation for the header file dolfin/mesh/Cell.h

.. _programmers_reference_cpp_mesh_cell:

Cell.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Cell

    *Parent class(es)*
    
        * :cpp:class:`MeshEntity`
        
    A Cell is a :cpp:class:`MeshEntity` of topological codimension 0.


    .. cpp:function:: Cell()
    
        Create empty cell


    .. cpp:function:: Cell(const Mesh& mesh, uint index)
    
        Create cell on given mesh with given index
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            index (uint)
                The index.


    .. cpp:function:: CellType::Type type() const
    
        Return type of cell


    .. cpp:function:: double orientation() const
    
        Compute orientation of cell
        
        *Returns*
            double
                Orientation of the cell (0 is right, 1 is left).


    .. cpp:function:: double volume() const
    
        Compute (generalized) volume of cell
        
        *Returns*
            double
                The volume of the cell.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Cell cell(mesh, 0);
                info("%g", cell.volume());
        
            output::
        
                0.5


    .. cpp:function:: double diameter() const
    
        Compute diameter of cell
        
        *Returns*
            double
                The diameter of the cell.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Cell cell(mesh, 0);
                info("%g", cell.diameter());
        
            output::
        
                1.41421


    .. cpp:function:: double normal(uint facet, uint i) const
    
        Compute component i of normal of given facet with respect to the cell
        
        *Arguments*
            facet (uint)
                Index of facet.
            i (uint)
                Component.
        
        *Returns*
            double
                Component i of the normal of the facet.


    .. cpp:function:: Point normal(uint facet) const
    
        Compute normal of given facet with respect to the cell
        
        *Arguments*
            facet (uint)
                Index of facet.
        
        *Returns*
            :cpp:class:`Point`
                Normal of the facet.


    .. cpp:function:: double facet_area(uint facet) const
    
        Compute the area/length of given facet with respect to the cell
        
        *Arguments*
            facet (uint)
                Index of the facet.
        
        *Returns*
            double
                Area/length of the facet.


    .. cpp:function:: void order(const MeshFunction<uint>* global_vertex_indices)
    
        Order entities locally
        
        *Arguments*
            global_vertex_indices (:cpp:class:`MeshFunction` <uint>)
                The global vertex indices.


    .. cpp:function:: bool ordered(const MeshFunction<uint>* global_vertex_indices) const
    
        Check if entities are ordered
        
        *Arguments*
            global_vertex_indices (:cpp:class:`MeshFunction` <uint>)
                The global vertex indices.
        
        *Returns*
            bool
                True if ordered.


.. cpp:class:: CellIterator

    *Parent class(es)*
    
        * :cpp:class:`MeshEntityIterator`
        
    A CellIterator is a MeshEntityIterator of topological codimension 0.


.. cpp:class:: CellFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A CellFunction is a MeshFunction of topological codimension 0.


