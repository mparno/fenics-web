
.. Documentation for the header file dolfin/mesh/MeshEditor.h

.. _programmers_reference_cpp_mesh_mesheditor:

MeshEditor.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEditor

    A simple mesh editor for creating simplicial meshes in 1D, 2D
    and 3D.


    .. cpp:function:: MeshEditor()
    
        Constructor


    .. cpp:function:: void open(Mesh& mesh, uint tdim, uint gdim)
    
        Open mesh of given topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            tdim (uint)
                The topological dimension.
            gdim (uint)
                The geometrical dimension.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
        


    .. cpp:function:: void open(Mesh& mesh, CellType::Type type, uint tdim, uint gdim)
    
        Open mesh of given cell type, topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            type (CellType::Type)
                Cell type.
            tdim (uint)
                The topological dimension.
            gdim (uint)
                The geometrical dimension.


    .. cpp:function:: void open(Mesh& mesh, std::string type, uint tdim, uint gdim)
    
        Open mesh of given cell type, topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            type (std::string)
                Cell type.
            tdim (uint)
                The topological dimension.
            gdim (uint)
                The geometrical dimension.


    .. cpp:function:: void init_vertices(uint num_vertices)
    
        Specify number of vertices
        
        *Arguments*
            num_vertices (uint)
                The number of vertices.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_vertices(4);
        


    .. cpp:function:: void init_higher_order_vertices(uint num_higher_order_vertices)
    
        Specify number of vertices
        
        *Arguments*
            num_higher_order_vertices (uint)
                The number of higher order vertices.


    .. cpp:function:: void init_cells(uint num_cells)
    
        Specify number of cells
        
        *Arguments*
            num_cells (uint)
                The number of cells.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_cells(2);
        


    .. cpp:function:: void init_higher_order_cells(uint num_higher_order_cells, uint num_higher_order_cell_dof)
    
        Specify number of cells
        
        *Arguments*
            num_higher_order_cells (uint)
                The number of higher order cells.
            num_higher_order_cell_dof (uint)
                The number of cell dofs.


    .. cpp:function:: void set_affine_cell_indicator(uint c, const std::string affine_str)
    
        Set boolean indicator inside MeshGeometry


    .. cpp:function:: void add_vertex(uint v, const Point& p)
    
        Add vertex v at given point p
        
        *Arguments*
            v (uint)
                The vertex (index).
            p (:cpp:class:`Point`)
                The point.


    .. cpp:function:: void add_vertex(uint v, double x)
    
        Add vertex v at given coordinate x
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.


    .. cpp:function:: void add_vertex(uint v, double x, double y)
    
        Add vertex v at given coordinate (x, y)
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.
        
        *Example*
            .. code-block:: c++
        
                MeshEditor editor;
                editor.add_vertex(0, 0.0, 0.0);
        


    .. cpp:function:: void add_vertex(uint v, double x, double y, double z)
    
        Add vertex v at given coordinate (x, y, z)
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.
            z (double)
                The z-coordinate.


    .. cpp:function:: void add_higher_order_vertex(uint v, const Point& p)
    
        Add vertex v at given point p
        
        *Arguments*
            v (uint)
                The vertex (index).
            p (:cpp:class:`Point`)
                The point.


    .. cpp:function:: void add_higher_order_vertex(uint v, double x)
    
        Add vertex v at given coordinate x
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.


    .. cpp:function:: void add_higher_order_vertex(uint v, double x, double y)
    
        Add vertex v at given coordinate (x, y)
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.


    .. cpp:function:: void add_higher_order_vertex(uint v, double x, double y, double z)
    
        Add vertex v at given coordinate (x, y, z)
        
        *Arguments*
            v (uint)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.
            z (double)
                The z-coordinate.


    .. cpp:function:: void add_cell(uint c, const std::vector<uint>& v)
    
        Add cell with given vertices
        
        *Arguments*
            c (uint)
                The cell (index).
            v (std::vector<uint>)
                The vertex indices


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1)
    
        Add cell (interval) with given vertices
        
        *Arguments*
            c (uint)
                The cell (index).
            v0 (uint)
                Index of the first vertex.
            v1 (uint)
                Index of the second vertex.


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1, uint v2)
    
        Add cell (triangle) with given vertices
        
        *Arguments*
            c (uint)
                The cell (index).
            v0 (uint)
                Index of the first vertex.
            v1 (uint)
                Index of the second vertex.
            v2 (uint)
                Index of the third vertex.
        
        *Example*
            .. code-block:: c++
        
                MeshEditor editor;
                editor.add_cell(0, 0, 1, 2);
        


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1, uint v2, uint v3)
    
        Add cell (tetrahedron) with given vertices
        
        *Arguments*
            c (uint)
                The cell (index).
            v0 (uint)
                Index of the first vertex.
            v1 (uint)
                Index of the second vertex.
            v2 (uint)
                Index of the third vertex.
            v3 (uint)
                Index of the fourth vertex.


    .. cpp:function:: void add_higher_order_cell_data(uint c, uint v0, uint v1, uint v2, uint v3, uint v4, uint v5)
    
        Add higher order cell data (assume P2 triangle for now)
        
        *Arguments*
            c (uint)
                The cell (index).
            v0 (uint)
                Index of the first vertex.
            v1 (uint)
                Index of the second vertex.
            v2 (uint)
                Index of the third vertex.
            v3 (uint)
                Index of the fourth vertex.
            v4 (uint)
                Index of the fifth vertex.
            v5 (uint)
                Index of the sixth vertex.


    .. cpp:function:: void close(bool order=true)
    
        Close mesh, finish editing, and order entities locally
        
        *Arguments*
            order (bool)
                Order entities locally if true. Default values is true.
        
        *Example*
            .. code-block:: c++
        
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                ...
                editor.close()
        


