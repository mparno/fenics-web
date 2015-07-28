
.. Documentation for the header file dolfin/mesh/MeshTopology.h

.. _programmers_reference_cpp_mesh_meshtopology:

MeshTopology.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshTopology

    MeshTopology stores the topology of a mesh, consisting of mesh
    entities and connectivity (incidence relations for the mesh
    entities). Note that the mesh entities don't need to be stored,
    only the number of entities and the connectivity. Any numbering
    scheme for the mesh entities is stored separately in a
    MeshFunction over the entities.
    
    A mesh entity e may be identified globally as a pair e = (dim,
    i), where dim is the topological dimension and i is the index of
    the entity within that topological dimension.


    .. cpp:function:: MeshTopology()
    
        Create empty mesh topology


    .. cpp:function:: MeshTopology(const MeshTopology& topology)
    
        Copy constructor


    .. cpp:function:: MeshTopology& operator= (const MeshTopology& topology)
    
        Assignment


    .. cpp:function:: std::size_t dim() const
    
        Return topological dimension


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return number of entities for given dimension


    .. cpp:function:: std::size_t size_global(std::size_t dim) const
    
        Return global number of entities for given dimension


    .. cpp:function:: std::size_t ghost_offset(std::size_t dim) const
    
        Return number of regular (non-ghost) entities
        or equivalently, the offset of where ghost entities begin


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: void clear(std::size_t d0, std::size_t d1)
    
        Clear data for given pair of topological dimensions


    .. cpp:function:: void init(std::size_t dim)
    
        Initialize topology of given maximum dimension


    .. cpp:function:: void init(std::size_t dim, std::size_t local_size, std::size_t global_size)
    
        Set number of local entities (local_size) and global entities
        (global_size) for given topological dimension dim


    .. cpp:function:: void init_global_indices(std::size_t dim, std::size_t size)
    
        Initialize storage for global entity numbering for entities of
        dimension dim


    .. cpp:function:: void init_ghost(std::size_t dim, std::size_t index)
    
        Initialise the offset index of ghost entities for this dimension


    .. cpp:function:: void set_global_index(std::size_t dim, std::size_t local_index, std::size_t global_index)
    
        Set global index for entity of dimension dim and with local
        index


    .. cpp:function:: const std::vector<std::size_t>& global_indices(std::size_t d) const
    
        Get local-to-global index map for entities of topological
        dimension d


    .. cpp:function:: bool have_global_indices(std::size_t dim) const
    
        Check if global indices are available for entities of
        dimension dim


    .. cpp:function:: bool have_shared_entities(unsigned int dim) const
    
        Check whether there are any shared entities calculated
        of dimension dim


    .. cpp:function:: std::map<unsigned int, std::set<unsigned int> >& shared_entities(unsigned int dim)
    
        Return map from shared entities (local index) to processes
        that share the entity


    .. cpp:function:: const std::map<unsigned int, std::set<unsigned int> >& shared_entities(unsigned int dim) const
    
        Return map from shared entities (local index) to process that
        share the entity (const version)


    .. cpp:function:: std::vector<unsigned int>& cell_owner()
    
        Return mapping from local ghost cell index to owning process
        Since ghost cells are at the end of the range, this is just
        a vector over those cells


    .. cpp:function:: const std::vector<unsigned int>& cell_owner() const
    
        Return mapping from local ghost cell index to owning process (const version)
        Since ghost cells are at the end of the range, this is just
        a vector over those cells


    .. cpp:function:: dolfin::MeshConnectivity& operator() (std::size_t d0, std::size_t d1)
    
        Return connectivity for given pair of topological dimensions


    .. cpp:function:: const dolfin::MeshConnectivity& operator() (std::size_t d0, std::size_t d1) const
    
        Return connectivity for given pair of topological dimensions


    .. cpp:function:: size_t hash() const
    
        Return hash based on the hash of cell-vertex connectivity


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


