
.. Documentation for the header file dolfin/mesh/MeshEditor.h

.. _programmers_reference_cpp_mesh_mesheditor:

MeshEditor.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEditor

    A simple mesh editor for creating simplicial meshes in 1D, 2D and 3D.


    .. cpp:function:: MeshEditor()
    
        Constructor


    .. cpp:function:: void open(Mesh& mesh, uint tdim, uint gdim)
    
        Open mesh of given topological and geometrical dimension


    .. cpp:function:: void open(Mesh& mesh, CellType::Type type, uint tdim, uint gdim)
    
        Open mesh of given cell type, topological and geometrical dimension


    .. cpp:function:: void open(Mesh& mesh, std::string type, uint tdim, uint gdim)
    
        Open mesh of given cell type, topological and geometrical dimension


    .. cpp:function:: void init_vertices(uint num_vertices)
    
        Specify number of vertices


    .. cpp:function:: void init_higher_order_vertices(uint num_higher_order_vertices)
    
        Specify number of vertices


    .. cpp:function:: void init_cells(uint num_cells)
    
        Specify number of cells


    .. cpp:function:: void init_higher_order_cells(uint num_higher_order_cells, uint num_higher_order_cell_dof)
    
        Specify number of cells


    .. cpp:function:: void set_affine_cell_indicator(uint c, const std::string affine_str)
    
        Set boolean indicator inside MeshGeometry


    .. cpp:function:: void add_vertex(uint v, const Point& p)
    
        Add vertex v at given point p


    .. cpp:function:: void add_vertex(uint v, double x)
    
        Add vertex v at given coordinate x


    .. cpp:function:: void add_vertex(uint v, double x, double y)
    
        Add vertex v at given coordinate (x, y)


    .. cpp:function:: void add_vertex(uint v, double x, double y, double z)
    
        Add vertex v at given coordinate (x, y, z)


    .. cpp:function:: void add_higher_order_vertex(uint v, const Point& p)
    
        Add vertex v at given point p


    .. cpp:function:: void add_higher_order_vertex(uint v, double x)
    
        Add vertex v at given coordinate x


    .. cpp:function:: void add_higher_order_vertex(uint v, double x, double y)
    
        Add vertex v at given coordinate (x, y)


    .. cpp:function:: void add_higher_order_vertex(uint v, double x, double y, double z)
    
        Add vertex v at given coordinate (x, y, z)


    .. cpp:function:: void add_cell(uint c, const std::vector<uint>& v)
    
        Add cell with given vertices


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1)
    
        Add cell (interval) with given vertices


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1, uint v2)
    
        Add cell (triangle) with given vertices


    .. cpp:function:: void add_cell(uint c, uint v0, uint v1, uint v2, uint v3)
    
        Add cell (tetrahedron) with given vertices


    .. cpp:function:: void add_higher_order_cell_data(uint c, uint v0, uint v1, uint v2, uint v3, uint v4, uint v5)
    
        Add higher order cell data (assume P2 triangle for now)


    .. cpp:function:: void close(bool order=true)
    
        Close mesh, finish editing, and order entities locally


