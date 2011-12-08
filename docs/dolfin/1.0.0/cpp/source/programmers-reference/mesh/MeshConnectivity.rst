
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


    .. cpp:function:: MeshConnectivity(uint d0, uint d1)
    
        Create empty connectivity between given dimensions (d0 -- d1)


    .. cpp:function:: MeshConnectivity(const MeshConnectivity& connectivity)
    
        Copy constructor


    .. cpp:function:: const MeshConnectivity& operator= (const MeshConnectivity& connectivity)
    
        Assignment


    .. cpp:function:: uint size() const
    
        Return total number of connections


    .. cpp:function:: uint size(uint entity) const
    
        Return number of connections for given entity


    .. cpp:function:: const uint* operator() (uint entity) const
    
        Return array of connections for given entity


    .. cpp:function:: const uint* operator() () const
    
        Return contiguous array of connections for all entities


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: void init(uint num_entities, uint num_connections)
    
        Initialize number of entities and number of connections (equal for all)


    .. cpp:function:: void init(std::vector<uint>& num_connections)
    
        Initialize number of entities and number of connections (individually)


    .. cpp:function:: void set(uint entity, uint connection, uint pos)
    
        Set given connection for given entity


    .. cpp:function:: void set(uint entity, const std::vector<uint>& connections)
    
        Set all connections for given entity


    .. cpp:function:: void set(uint entity, uint* connections)
    
        Set all connections for given entity


    .. cpp:function:: void set(const std::vector<std::vector<uint> >& connectivity)
    
        Set all connections for all entities


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


