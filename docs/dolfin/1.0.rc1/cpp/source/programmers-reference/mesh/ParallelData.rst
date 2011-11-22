
.. Documentation for the header file dolfin/mesh/ParallelData.h

.. _programmers_reference_cpp_mesh_paralleldata:

ParallelData.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ParallelData

    This class stores auxiliary mesh data for parallel computing.


    .. cpp:function:: ParallelData(const Mesh& mesh)
    
        Constructor


    .. cpp:function:: ParallelData(const ParallelData& data)
    
        Copy constructor


    .. cpp:function:: bool have_global_entity_indices(uint d) const
    
        Return true if global indices have been computed for entity of
        dimension d


    .. cpp:function:: MeshFunction<unsigned int>& global_entity_indices(uint d)
    
        Return global indices (local-to-global) for entity of dimension d


    .. cpp:function:: const MeshFunction<unsigned int>& global_entity_indices(uint d) const
    
        Return global indices (local-to-global) for entity of dimension d (const version)


    .. cpp:function:: std::vector<uint> global_entity_indices_as_vector(uint d) const
    
        Return global indices (local-to-global) for entity of dimension d in a vector


    .. cpp:function:: const std::map<unsigned int, unsigned int>& global_to_local_entity_indices(uint d)
    
        Return global-to-local indices for entity of dimension d


    .. cpp:function:: const std::map<unsigned int, unsigned int>& global_to_local_entity_indices(uint d) const
    
        Return global-to-local indices for entity of dimension d (const version)


    .. cpp:function:: std::map<unsigned int, std::vector<unsigned int> >& shared_vertices()
    
        FIXME: Add description and use better name


    .. cpp:function:: const std::map<unsigned int, std::vector<unsigned int> >& shared_vertices() const
    
        FIXME: Add description and use better name


    .. cpp:function:: MeshFunction<bool>& exterior_facet()
    
        Return MeshFunction that is true for globally exterior facets,
        false otherwise


    .. cpp:function:: const MeshFunction<bool>& exterior_facet() const
    
        Return MeshFunction that is true for globally exterior facets,
        false otherwise (const version)


