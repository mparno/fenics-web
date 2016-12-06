
.. Documentation for the header file dolfin/graph/SCOTCH.h

.. _programmers_reference_cpp_graph_scotch:

SCOTCH.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SCOTCH

    This class provides an interface to SCOTCH-PT (parallel version)


    .. cpp:function:: static void compute_partition( const MPI_Comm mpi_comm, std::vector<int>& cell_partition, std::map<std::int64_t, std::vector<int>>& ghost_procs, const boost::multi_array<std::int64_t, 2>& cell_vertices, const std::vector<std::size_t>& cell_weight, const std::int64_t num_global_vertices, const std::int64_t num_global_cells, const CellType& cell_type)
    
        Compute cell partition from local mesh data.  The vector
        cell_partition contains the desired destination process
        numbers for each cell.  Cells shared on multiple processes
        have an entry in ghost_procs pointing to the set of sharing
        process numbers.


    .. cpp:function:: static std::vector<int> compute_gps(const Graph& graph, std::size_t num_passes=5)
    
        Compute reordering (map[old] -> new) using
        Gibbs-Poole-Stockmeyer (GPS) re-ordering


