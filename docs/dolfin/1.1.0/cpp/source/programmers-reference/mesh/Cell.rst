
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


    .. cpp:function:: Cell(const Mesh& mesh, std::size_t index)
    
        Create cell on given mesh with given index
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            index (std::size_t)
                The index.


    .. cpp:function:: CellType::Type type() const
    
        Return type of cell


    .. cpp:function:: std::size_t orientation() const
    
        Compute orientation of cell
        
        *Returns*
            std::size_t
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


    .. cpp:function:: double normal(std::size_t facet, std::size_t i) const
    
        Compute component i of normal of given facet with respect to the cell
        
        *Arguments*
            facet (std::size_t)
                Index of facet.
            i (std::size_t)
                Component.
        
        *Returns*
            double
                Component i of the normal of the facet.


    .. cpp:function:: Point normal(std::size_t facet) const
    
        Compute normal of given facet with respect to the cell
        
        *Arguments*
            facet (std::size_t)
                Index of facet.
        
        *Returns*
            :cpp:class:`Point`
                Normal of the facet.


    .. cpp:function:: double facet_area(std::size_t facet) const
    
        Compute the area/length of given facet with respect to the cell
        
        *Arguments*
            facet (std::size_t)
                Index of the facet.
        
        *Returns*
            double
                Area/length of the facet.


    .. cpp:function:: void order(const std::vector<std::size_t>& local_to_global_vertex_indices)
    
        Order entities locally
        
        *Arguments*
            global_vertex_indices (:cpp:class:`MeshFunction` <std::size_t>)
                The global vertex indices.


    .. cpp:function:: bool ordered(const std::vector<std::size_t>& local_to_global_vertex_indices) const
    
        Check if entities are ordered
        
        *Arguments*
            global_vertex_indices (:cpp:class:`MeshFunction` <std::size_t>)
                The global vertex indices.
        
        *Returns*
            bool
                True if ordered.


.. cpp:class:: CellFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A CellFunction is a MeshFunction of topological codimension 0.


