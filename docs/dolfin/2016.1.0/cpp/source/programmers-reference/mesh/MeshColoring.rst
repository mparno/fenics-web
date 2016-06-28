
.. Documentation for the header file dolfin/mesh/MeshColoring.h

.. _programmers_reference_cpp_mesh_meshcoloring:

MeshColoring.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshColoring

    This class computes colorings for a local mesh. It supports
    vertex, edge, and facet-based colorings.


    .. cpp:function:: static const std::vector<std::size_t>& color_cells(Mesh& mesh, std::string coloring_type)
    
        Color the cells of a mesh for given coloring type, which can
        be one of "vertex", "edge" or "facet". Coloring is saved in
        the mesh topology


    .. cpp:function:: static const std::vector<std::size_t>& color(Mesh& mesh, const std::vector<std::size_t>& coloring_type)
    
        Color the cells of a mesh for given coloring type specified by
        topological dimension, which can be one of 0, 1 or D -
        1. Coloring is saved in the mesh topology


    .. cpp:function:: static std::size_t compute_colors(const Mesh& mesh, std::vector<std::size_t>& colors, const std::vector<std::size_t>& coloring_type)
    
        Compute cell colors for given coloring type specified by
        topological dimension, which can be one of 0, 1 or D - 1.


    .. cpp:function:: static CellFunction<std::size_t> cell_colors(std::shared_ptr<const Mesh> mesh, std::string coloring_type)
    
        Return a MeshFunction with the cell colors (used for
        visualisation)


    .. cpp:function:: static CellFunction<std::size_t> cell_colors(std::shared_ptr<const Mesh> mesh, std::vector<std::size_t> coloring_type)
    
        Return a MeshFunction with the cell colors (used for visualisation)


    .. cpp:function:: static std::size_t type_to_dim(std::string coloring_type, const Mesh& mesh)
    
        Convert coloring type to topological dimension


