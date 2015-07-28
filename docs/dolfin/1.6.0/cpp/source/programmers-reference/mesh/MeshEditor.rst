
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


    .. cpp:function:: void open(Mesh& mesh, std::size_t tdim, std::size_t gdim)
    
        Open mesh of given topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            tdim (std::size_t)
                The topological dimension.
            gdim (std::size_t)
                The geometrical dimension.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
        


    .. cpp:function:: void open(Mesh& mesh, CellType::Type type, std::size_t tdim, std::size_t gdim)
    
        Open mesh of given cell type, topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            type (CellType::Type)
                Cell type.
            tdim (std::size_t)
                The topological dimension.
            gdim (std::size_t)
                The geometrical dimension.


    .. cpp:function:: void open(Mesh& mesh, std::string type, std::size_t tdim, std::size_t gdim)
    
        Open mesh of given cell type, topological and geometrical dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to open.
            type (std::string)
                Cell type.
            tdim (std::size_t)
                The topological dimension.
            gdim (std::size_t)
                The geometrical dimension.


    .. cpp:function:: void init_vertices(std::size_t num_vertices)
    
        Specify number of vertices (serial version)
        
        *Arguments*
            num_vertices (std::size_t)
                The number of vertices.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_vertices(9);
        


    .. cpp:function:: void init_vertices_global(std::size_t num_local_vertices, std::size_t num_global_vertices)
    
        Specify number of vertices (distributed version)
        
        *Arguments*
            num_local_vertices (std::size_t)
                The number of vertices on this process.
            num_global_vertices (std::size_t)
                The number of vertices in distributed mesh.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_vertices(4, 8);
        


    .. cpp:function:: void init_cells(std::size_t num_cells)
    
        Specify number of cells (serial version)
        
        *Arguments*
            num_cells (std::size_t)
                The number of cells.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_cells(8);


    .. cpp:function:: void init_cells_global(std::size_t num_local_cells, std::size_t num_global_cells)
    
        Specify number of cells (distributed version)
        
        *Arguments*
            num_local_cells (std::size_t)
                The number of local cells.
            num_global_cells (std::size_t)
                The number of cells in distributed mesh.
        
        *Example*
            .. code-block:: c++
        
                Mesh mesh;
                MeshEditor editor;
                editor.open(mesh, 2, 2);
                editor.init_cells(2, 6);
        


    .. cpp:function:: void add_vertex(std::size_t index, const Point& p)
    
        Add vertex v at given point p
        
        *Arguments*
            index (std::size_t)
                The vertex (index).
            p (:cpp:class:`Point`)
                The point.


    .. cpp:function:: void add_vertex(std::size_t index, const std::vector<double>& x)
    
        Add vertex v at given coordinate x
        
        *Arguments*
            index (std::size_t)
                The vertex (index).
            x (std::vector<double>)
                The x-coordinates.


    .. cpp:function:: void add_vertex(std::size_t index, double x)
    
        Add vertex v at given point x (for a 1D mesh)
        
        *Arguments*
            index (std::size_t)
                The vertex (index).
            x (double)
                The x-coordinate.


    .. cpp:function:: void add_vertex(std::size_t index, double x, double y)
    
        Add vertex v at given point (x, y) (for a 2D mesh)
        
        *Arguments*
            index (std::size_t)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.


    .. cpp:function:: void add_vertex(std::size_t index, double x, double y, double z)
    
        Add vertex v at given point (x, y, z) (for a 3D mesh)
        
        *Arguments*
            index (std::size_t)
                The vertex (index).
            x (double)
                The x-coordinate.
            y (double)
                The y-coordinate.
            z (double)
                The z-coordinate.


    .. cpp:function:: void add_vertex_global(std::size_t local_index, std::size_t global_index, const Point& p)
    
        Add vertex v at given point p
        
        *Arguments*
            local_index (std::size_t)
                The vertex (local index).
            global_index (std::size_t)
                The vertex (global_index).
            p (:cpp:class:`Point`)
                The point.


    .. cpp:function:: void add_vertex_global(std::size_t local_index, std::size_t global_index, const std::vector<double>& x)
    
        Add vertex v at given coordinate x
        
        *Arguments*
            local_index (std::size_t)
                The vertex (local index).
            global_index (std::size_t)
                The vertex (global_index).
            x (std::vector<double>)
                The x-coordinates.


    .. cpp:function:: void add_cell(std::size_t c, std::size_t v0, std::size_t v1)
    
        Add cell with given vertices (1D)
        
        *Arguments*
            c (std::size_t)
                The cell (index).
            v0 (std::vector<std::size_t>)
                The first vertex (local index).
            v1 (std::vector<std::size_t>)
                The second vertex (local index).


    .. cpp:function:: void add_cell(std::size_t c, std::size_t v0, std::size_t v1, std::size_t v2)
    
        Add cell with given vertices (2D)
        
        *Arguments*
            c (std::size_t)
                The cell (index).
            v0 (std::vector<std::size_t>)
                The first vertex (local index).
            v1 (std::vector<std::size_t>)
                The second vertex (local index).
            v2 (std::vector<std::size_t>)
                The third vertex (local index).


    .. cpp:function:: void add_cell(std::size_t c, std::size_t v0, std::size_t v1, std::size_t v2, std::size_t v3)
    
        Add cell with given vertices (3D)
        
        *Arguments*
            c (std::size_t)
                The cell (index).
            v0 (std::vector<std::size_t>)
                The first vertex (local index).
            v1 (std::vector<std::size_t>)
                The second vertex (local index).
            v2 (std::vector<std::size_t>)
                The third vertex (local index).
            v3 (std::vector<std::size_t>)
                The fourth vertex (local index).


    .. cpp:function:: void add_cell(std::size_t c, const std::vector<std::size_t>& v)
    
        Add cell with given vertices (non-templated version for Python
        interface)
        
        *Arguments*
            c (std::size_t)
                The cell (index).
            v (std::vector<std::size_t>)
                The vertex indices (local indices)


    .. cpp:function:: void add_cell(std::size_t c, const T& v)
    
        Add cell with given vertices
        
        *Arguments*
            c (std::size_t)
                The cell (index).
            v (typename T)
                The vertex indices (local indices)


    .. cpp:function:: void add_cell(std::size_t local_index, std::size_t global_index, const T& v)
    
        Add cell with given vertices
        
        *Arguments*
            local_index (std::size_t)
                The cell (index).
            global_index (std::size_t)
                The global (user) cell index.
            v (std::vector<std::size_t>)
                The vertex indices (local indices)


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
        


