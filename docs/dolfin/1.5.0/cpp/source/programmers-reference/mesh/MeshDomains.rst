
.. Documentation for the header file dolfin/mesh/MeshDomains.h

.. _programmers_reference_cpp_mesh_meshdomains:

MeshDomains.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshDomains

    The class :cpp:class:`MeshDomains` stores the division of a :cpp:class:`Mesh` into
    subdomains. For each topological dimension 0 <= d <= D, where D
    is the topological dimension of the :cpp:class:`Mesh`, a set of integer
    markers are stored for a subset of the entities of dimension d,
    indicating for each entity in the subset the number of the
    subdomain. It should be noted that the subset does not need to
    contain all entities of any given dimension; entities not
    contained in the subset are "unmarked".


    .. cpp:function:: MeshDomains()
    
        Create empty mesh domains


    .. cpp:function:: std::size_t max_dim() const
    
        Return maximum topological dimension of stored markers


    .. cpp:function:: std::size_t num_marked(std::size_t dim) const
    
        Return number of marked entities of given dimension


    .. cpp:function:: bool is_empty() const
    
        Check whether domain data is empty


    .. cpp:function:: std::map<std::size_t, std::size_t>& markers(std::size_t dim)
    
        Get subdomain markers for given dimension (shared pointer
        version)


    .. cpp:function:: const std::map<std::size_t, std::size_t>& markers(std::size_t dim) const
    
        Get subdomain markers for given dimension (const shared
        pointer version)


    .. cpp:function:: bool set_marker(std::pair<std::size_t, std::size_t> marker, std::size_t dim)
    
        Set marker (entity index, marker value) of a given dimension
        d. Returns true if a new key is inserted, false otherwise.


    .. cpp:function:: std::size_t get_marker(std::size_t entity_index, std::size_t dim) const
    
        Get marker (entity index, marker value) of a given dimension
        d. Throws an error if marker does not exist.


    .. cpp:function:: const MeshDomains& operator= (const MeshDomains& domains)
    
        Assignment operator


    .. cpp:function:: void init(std::size_t dim)
    
        Initialize mesh domains for given topological dimension


    .. cpp:function:: void clear()
    
        Clear all data


