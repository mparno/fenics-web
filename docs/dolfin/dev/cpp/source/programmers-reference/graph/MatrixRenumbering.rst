
.. Documentation for the header file dolfin/graph/MatrixRenumbering.h

.. _programmers_reference_cpp_graph_matrixrenumbering:

MatrixRenumbering.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MatrixRenumbering

    This class computes re-ordering based on a SparsityPattern graph
    representation of a sparse matrix. It uses Zoltan, which is part of
    Trilinos.


    .. cpp:function:: int num_global_objects() const
    
        Number of global graph vertices


    .. cpp:function:: int num_local_objects() const
    
        Number of local graph vertices


    .. cpp:function:: void num_edges_per_vertex(std::vector<uint>& num_edges) const
    
        Number of edges per vertex


    .. cpp:function:: const std::vector<std::vector<uint> > edges() const
    
        Vertex edges


