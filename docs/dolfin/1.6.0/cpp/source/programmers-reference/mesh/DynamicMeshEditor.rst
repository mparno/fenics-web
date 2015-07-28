
.. Documentation for the header file dolfin/mesh/DynamicMeshEditor.h

.. _programmers_reference_cpp_mesh_dynamicmesheditor:

DynamicMeshEditor.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DynamicMeshEditor

    This class provides an interface for dynamic editing of meshes,
    that is, when the number of vertices and cells are not known
    a priori. If the number of vertices and cells are known a priori,
    it is more efficient to use the default editor MeshEditor.


    .. cpp:function:: DynamicMeshEditor()
    
        Constructor


    .. cpp:function:: void open(Mesh& mesh, CellType::Type type, std::size_t tdim, std::size_t gdim, std::size_t num_global_vertices, std::size_t num_global_cells)
    
        Open mesh of given cell type, topological and geometrical dimension


    .. cpp:function:: void open(Mesh& mesh, std::string type, std::size_t tdim, std::size_t gdim, std::size_t num_global_vertices, std::size_t num_global_cells)
    
        Open mesh of given cell type, topological and geometrical dimension


    .. cpp:function:: void add_vertex(std::size_t v, const Point& p)
    
        Add vertex v at given point p


    .. cpp:function:: void add_vertex(std::size_t v, double x)
    
        Add vertex v at given coordinate x


    .. cpp:function:: void add_vertex(std::size_t v, double x, double y)
    
        Add vertex v at given coordinate (x, y)


    .. cpp:function:: void add_vertex(std::size_t v, double x, double y, double z)
    
        Add vertex v at given coordinate (x, y, z)


    .. cpp:function:: void add_cell(std::size_t c, const std::vector<std::size_t>& v)
    
        Add cell with given vertices


    .. cpp:function:: void add_cell(std::size_t c, std::size_t v0, std::size_t v1)
    
        Add cell (interval) with given vertices


    .. cpp:function:: void add_cell(std::size_t c,  std::size_t v0, std::size_t v1, std::size_t v2)
    
        Add cell (triangle) with given vertices


    .. cpp:function:: void add_cell(std::size_t c, std::size_t v0, std::size_t v1, std::size_t v2, std::size_t v3)
    
        Add cell (tetrahedron) with given vertices


    .. cpp:function:: void close(bool order=false)
    
        Close mesh, finish editing, and order entities locally


