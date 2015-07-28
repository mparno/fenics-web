
.. Documentation for the header file dolfin/graph/BoostGraphOrdering.h

.. _programmers_reference_cpp_graph_boostgraphordering:

BoostGraphOrdering.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BoostGraphOrdering

    This class computes graph re-orderings. It uses Boost Graph.


    .. cpp:function:: static std::vector<int> compute_cuthill_mckee(const Graph& graph, bool reverse=false)
    
        Compute re-ordering (map[old] -> new) using Cuthill-McKee algorithm


    .. cpp:function:: static std::vector<int> compute_cuthill_mckee(const std::set<std::pair<std::size_t, std::size_t> >& edges, std::size_t size, bool reverse=false)
    
        Compute re-ordering (map[old] -> new) using Cuthill-McKee algorithm


