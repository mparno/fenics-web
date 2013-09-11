
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
                Orientation of the cell (0 is 'up'/'right', 1 is 'down'/'left')


    .. cpp:function:: std::size_t orientation(const Point& up) const
    
        Compute orientation of cell relative to given 'up' direction
        
        *Arguments*
            up (:cpp:class:`Point`)
                The direction defined as 'up'
        
        *Returns*
            std::size_t
                Orientation of the cell (0 is 'same', 1 is 'opposite')


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


    .. cpp:function:: double inradius() const
    
        Compute inradius of cell
        
        *Returns*
            double
                Radius of the sphere inscribed in the cell.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Cell cell(mesh, 0);
                info("%g", cell.inradius());
        
            output::
        
                0.29289


    .. cpp:function:: double radius_ratio() const
    
        Compute ratio of inradius to circumradius times dim for cell.
        Useful as cell quality measure. Returns 1. for equilateral
        and 0. for degenerate cell.
        See Jonathan Richard Shewchuk: What Is a Good Linear Finite Element?,
        online: http://www.cs.berkeley.edu/~jrs/papers/elemj.pdf
        
        *Returns*
            double
                cell_dimension * inradius / circumradius
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Cell cell(mesh, 0);
                info("%g", cell.radius_ratio());
        
            output::
        
                0.828427


    .. cpp:function:: double squared_distance(const Point& point)
    
        Compute squared distance to given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
        *Returns*
            double
                The squared distance to the point.


    .. cpp:function:: double distance(const Point& point)
    
        Compute distance to given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
        *Returns*
            double
                The distance to the point.


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


    .. cpp:function:: Point cell_normal() const
    
        Compute normal to cell itself (viewed as embedded in 3D)
        
        *Returns*
            :cpp:class:`Point`
                Normal of the cell


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
                True iff ordered.


    .. cpp:function:: bool contains(const Point& point) const
    
        Check whether given point is contained in cell
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point to be checked.
        
        *Returns*
            bool
                True iff point is contained in cell.


.. cpp:class:: CellFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A CellFunction is a MeshFunction of topological codimension 0.


