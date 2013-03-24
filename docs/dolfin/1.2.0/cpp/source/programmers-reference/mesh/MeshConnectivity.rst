
.. Documentation for the header file dolfin/mesh/MeshConnectivity.h

.. _programmers_reference_cpp_mesh_meshconnectivity:

MeshConnectivity.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshConnectivity

    Mesh connectivity stores a sparse data structure of connections
    (incidence relations) between mesh entities for a fixed pair of
    topological dimensions.
    
    The connectivity can be specified either by first giving the
    number of entities and the number of connections for each entity,
    which may either be equal for all entities or different, or by
    giving the entire (sparse) connectivity pattern.


    .. cpp:function:: MeshConnectivity(std::size_t d0, std::size_t d1)
    
        Create empty connectivity between given dimensions (d0 -- d1)


    .. cpp:function:: MeshConnectivity(const MeshConnectivity& connectivity)
    
        Copy constructor


    .. cpp:function:: const MeshConnectivity& operator= (const MeshConnectivity& connectivity)
    
        Assignment


    .. cpp:function:: bool empty() const
    
        Return true if the total number of connections is equal to zero


    .. cpp:function:: std::size_t size() const
    
        Return total number of connections


    .. cpp:function:: std::size_t size(std::size_t entity) const
    
        Return number of connections for given entity


    .. cpp:function:: std::size_t size_global(std::size_t entity) const
    
        Return global number of connections for given entity


    .. cpp:function:: const unsigned int* operator() (std::size_t entity) const
    
        Return array of connections for given entity


    .. cpp:function:: const std::vector<unsigned int>& operator() () const
    
        Return contiguous array of connections for all entities


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: void init(std::size_t num_entities, std::size_t num_connections)
    
        Initialize number of entities and number of connections (equal
        for all)


    .. cpp:function:: void init(std::vector<std::size_t>& num_connections)
    
        Initialize number of entities and number of connections
        (individually)


    .. cpp:function:: void set(std::size_t entity, std::size_t connection, std::size_t pos)
    
        Set given connection for given entity


    .. cpp:function:: void set(std::size_t entity, const std::vector<std::size_t>& connections)
    
        Set all connections for given entity


    .. cpp:function:: void set(std::size_t entity, std::size_t* connections)
    
        Set all connections for given entity


    .. cpp:function:: void set(const std::vector<T>& connections)
    
        Set all connections for all entities (T is a container, e.g.
        a std::vector<std::size_t>, std::set<std::size_t>, etc)


    .. cpp:function:: void set_global_size(const std::vector<unsigned int>& num_global_connections)
    
        Set global number of connections for all local entities


    .. cpp:function:: std::size_t hash() const
    
        Hash of connections


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


