
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


    .. cpp:function:: static const MeshFunction<uint>& color_cells(Mesh& mesh, std::string coloring_type)
    
        Color the cells of a mesh for given coloring type, which can
        be one of "vertex", "edge" or "facet".


    .. cpp:function:: static const MeshFunction<uint>& color(Mesh& mesh, std::vector<uint> coloring_type)
    
        Color the cells of a mesh for given coloring type specified by
        topological dimension, which can be one of 0, 1 or D - 1.


    .. cpp:function:: static uint compute_colors(MeshFunction<uint>& colors, const std::vector<uint> coloring_type)
    
        Compute cell colors for given coloring type specified by
        topological dimension, which can be one of 0, 1 or D - 1.


    .. cpp:function:: static uint type_to_dim(std::string coloring_type, const Mesh& mesh)
    
        Convert coloring type to topological dimension


