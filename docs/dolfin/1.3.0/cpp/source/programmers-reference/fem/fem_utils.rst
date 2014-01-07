
.. Documentation for the header file dolfin/fem/fem_utils.h

.. _programmers_reference_cpp_fem_fem_utils:

fem_utils.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: std::vector<std::size_t> dof_to_vertex_map(const FunctionSpace& space)

    Return a map between dofs indices and vertex indices
    
    Only works for FunctionSpace with dofs exclusively on vertices.
    For MixedFunctionSpaces vertex index is offset with the number
    of dofs per vertex. In parallel the returned map only maps local
    (to processor) dofs.
    
    *Arguments*
        space (:cpp:class:`FunctionSpace`)
            The FunctionSpace for what the dof to vertex map should be computed for
    
    *Returns*
        std::vector<std::size_t>
            The dof to vertex map


.. cpp:function:: std::vector<dolfin::la_index> vertex_to_dof_map(const FunctionSpace& space)

    Return a map between vertex indices and dofs indices
    
    Only works for FunctionSpace with dofs exclusively on vertices.
    For MixedFunctionSpaces dof index is offset with the number of
    dofs per vertex.
    
    *Arguments*
        space (:cpp:class:`FunctionSpace`)
            The FunctionSpace for what the vertex to dof map should be computed for
    
    *Returns*
        std::vector<dolfin::la_index>
            The vertex to dof map


