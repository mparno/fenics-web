
.. Documentation for the header file dolfin/graph/GraphBuilder.h

.. _programmers_reference_cpp_graph_graphbuilder:

GraphBuilder.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GraphBuilder

    This class builds a Graph corresponding to various objects


    .. cpp:function:: static Graph local_graph(const Mesh& mesh, const GenericDofMap& dofmap0, const GenericDofMap& dofmap1)
    
        Build local graph from dofmap


    .. cpp:function:: static Graph local_graph(const Mesh& mesh, const std::vector<std::size_t>& coloring_type)
    
        Build local graph from mesh (general version)


    .. cpp:function:: static Graph local_graph(const Mesh& mesh, std::size_t dim0, std::size_t dim1)
    
        Build local graph (specialized version)


    .. cpp:function:: static std::pair<std::int32_t, std::int32_t> compute_dual_graph(const MPI_Comm mpi_comm, const boost::multi_array<std::int64_t, 2>& cell_vertices, const CellType& cell_type, const std::int64_t num_global_vertices, std::vector<std::vector<std::size_t>>& local_graph, std::set<std::int64_t>& ghost_vertices)
    
        Build distributed dual graph (cell-cell connections) from
        minimal mesh data, and return (num local edges, num
        non-local edges)


