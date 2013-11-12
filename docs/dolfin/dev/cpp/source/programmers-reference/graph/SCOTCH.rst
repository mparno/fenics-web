
.. Documentation for the header file dolfin/graph/SCOTCH.h

.. _programmers_reference_cpp_graph_scotch:

SCOTCH.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SCOTCH

    This class proivdes an interface to SCOTCH-PT (parallel version)


    .. cpp:function:: static std::vector<std::size_t> compute_gps(const Graph& graph, std::size_t num_passes=5)
    
        Compute reordering (map[old] -> new) using
        Gibbs-Poole-Stockmeyer re-ordering


