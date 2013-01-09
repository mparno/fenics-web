
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


    .. cpp:function:: static BoostBidirectionalGraph local_boost_graph(const Mesh& mesh, const std::vector<std::size_t>& coloring_type)
    
        Build local Boost graph (general version)


    .. cpp:function:: static Graph local_graph(const Mesh& mesh, std::size_t dim0, std::size_t dim1)
    
        Build local graph (specialized version)


    .. cpp:function:: static BoostBidirectionalGraph local_boost_graph(const Mesh& mesh, std::size_t dim0, std::size_t dim1)
    
        Build local Boost graph (specialized version)


    .. cpp:function:: static void compute_dual_graph(const LocalMeshData& mesh_data, std::vector<std::set<std::size_t> >& local_graph, std::set<std::size_t>& ghost_vertices)
    
        Build distributed dual graph for mesh. This function is slower
        than 'compute_dual_graph_small' for low processes counts, but
        scales much better with increasing process count.


    .. cpp:function:: static void compute_dual_graph_small(const LocalMeshData& mesh_data, std::vector<std::set<std::size_t> >& local_graph, std::set<std::size_t>& ghost_vertices)
    
        Build distributed dual graph for mesh. This function is very
        fast for a small number of processes, but does not scale for
        increasing process count.


