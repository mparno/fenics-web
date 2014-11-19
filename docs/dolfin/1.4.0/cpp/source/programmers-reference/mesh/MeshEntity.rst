
.. Documentation for the header file dolfin/mesh/MeshEntity.h

.. _programmers_reference_cpp_mesh_meshentity:

MeshEntity.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEntity

    A MeshEntity represents a mesh entity associated with
    a specific topological dimension of some :cpp:class:`Mesh`.


    .. cpp:function:: MeshEntity()
    
        Default Constructor


    .. cpp:function:: MeshEntity(const Mesh& mesh, std::size_t dim, std::size_t index)
    
        Constructor
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            dim (std::size_t)
                The topological dimension.
            index (std::size_t)
                The index.


    .. cpp:function:: void init(const Mesh& mesh, std::size_t dim, std::size_t index)
    
        Initialize mesh entity with given data
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            dim (std::size_t)
                The topological dimension.
            index (std::size_t)
                The index.


    .. cpp:function:: bool operator==(const MeshEntity& e) const
    
        Comparision Operator
        
        *Arguments*
            another (:cpp:class:`MeshEntity`)
                Another mesh entity
        
        *Returns*
            bool
                True if the two mesh entities are equal.


    .. cpp:function:: bool operator!=(const MeshEntity& e) const
    
        Comparision Operator
        
        *Arguments*
            another (MeshEntity)
                Another mesh entity.
        
        *Returns*
            bool
                True if the two mesh entities are NOT equal.


    .. cpp:function:: const Mesh& mesh() const
    
        Return mesh associated with mesh entity
        
        *Returns*
            :cpp:class:`Mesh`
                The mesh.


    .. cpp:function:: std::size_t dim() const
    
        Return topological dimension
        
        *Returns*
            std::size_t
                The dimension.


    .. cpp:function:: std::size_t index() const
    
        Return index of mesh entity
        
        *Returns*
            std::size_t
                The index.


    .. cpp:function:: std::size_t global_index() const
    
        Return global index of mesh entity
        
        *Returns*
            std::size_t
                The global index. Set to
                std::numerical_limits<std::size_t>::max() if global index
                has not been computed


    .. cpp:function:: std::size_t num_entities(std::size_t dim) const
    
        Return local number of incident mesh entities of given
        topological dimension
        
        *Arguments*
            dim (std::size_t)
                The topological dimension.
        
        *Returns*
            std::size_t
        The number of local incident MeshEntity objects of given
        dimension.


    .. cpp:function:: std::size_t num_global_entities(std::size_t dim) const
    
        Return global number of incident mesh entities of given
        topological dimension
        
        *Arguments*
            dim (std::size_t)
                The topological dimension.
        
        *Returns*
            std::size_t
                The number of global incident MeshEntity objects of given
                dimension.


    .. cpp:function:: const unsigned int* entities(std::size_t dim) const
    
        Return array of indices for incident mesh entitites of given
        topological dimension
        
        *Arguments*
            dim (std::size_t)
                The topological dimension.
        
        *Returns*
            std::size_t
                The index for incident mesh entities of given dimension.


    .. cpp:function:: std::size_t mesh_id() const
    
        Return unique mesh ID
        
        *Returns*
            std::size_t
                The unique mesh ID.


    .. cpp:function:: bool incident(const MeshEntity& entity) const
    
        Check if given entity is incident
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The entity.
        
        *Returns*
            bool
                True if the given entity is incident


    .. cpp:function:: std::size_t index(const MeshEntity& entity) const
    
        Compute local index of given incident entity (error if not
        found)
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The mesh entity.
        
        *Returns*
            std::size_t
                The local index of given entity.


    .. cpp:function:: Point midpoint() const
    
        Compute midpoint of cell
        
        *Returns*
            :cpp:class:`Point`
                The midpoint of the cell.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


