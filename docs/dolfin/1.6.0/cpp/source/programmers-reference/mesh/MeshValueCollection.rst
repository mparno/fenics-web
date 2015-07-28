
.. Documentation for the header file dolfin/mesh/MeshValueCollection.h

.. _programmers_reference_cpp_mesh_meshvaluecollection:

MeshValueCollection.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshValueCollection

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    The MeshValueCollection class can be used to store data
    associated with a subset of the entities of a mesh of a given
    topological dimension. It differs from the MeshFunction class in
    two ways. First, data does not need to be associated with all
    entities (only a subset). Second, data is associated with
    entities through the corresponding cell index and local entity
    number (relative to the cell), not by global entity index, which
    means that data may be stored robustly to file.


    .. cpp:function:: MeshValueCollection()
    
        Create empty mesh value collection
        


    .. cpp:function:: explicit MeshValueCollection(std::shared_ptr<const Mesh> mesh)
    
        Create an empty mesh value collection on a given mesh
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: explicit MeshValueCollection(const MeshFunction<T>& mesh_function)
    
        Create a mesh value collection from a MeshFunction
        
        *Arguments*
            mesh_function (:cpp:class:`MeshFunction` <T>)
                The mesh function for creating a MeshValueCollection.


    .. cpp:function:: MeshValueCollection(const Mesh& mesh, std::size_t dim)
    
        Create a mesh value collection of entities of given dimension
        on a given mesh
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh associated with the collection.
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: MeshValueCollection(std::shared_ptr<const Mesh> mesh, std::size_t dim)
    
        Create a mesh value collection of entities of given dimension
        on a given mesh (shared_ptr version)
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh associated with the collection.
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: MeshValueCollection(const Mesh& mesh, const std::string filename)
    
        Create a mesh value collection from a file.
        
        *Arguments*
            mesh (Mesh)
                A mesh associated with the collection. The mesh is used to
                map collection values to the appropriate process.
            filename (std::string)
                The XML file name.
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: MeshValueCollection<T>& operator=(const MeshFunction<T>& mesh_function)
    
        Assignment operator
        
        *Arguments*
            mesh_function (:cpp:class:`MeshFunction`)
                A :cpp:class:`MeshFunction` object used to construct a
                MeshValueCollection.


    .. cpp:function:: MeshValueCollection<T>& operator=(const MeshValueCollection<T>& mesh_value_collection)
    
        Assignment operator
        
        *Arguments*
            mesh_value_collection (:cpp:class:`MeshValueCollection`)
                A :cpp:class:`MeshValueCollection` object used to construct a
                MeshValueCollection.


    .. cpp:function:: void init(const Mesh& mesh, std::size_t dim)
    
        Initialise MeshValueCollection with mesh and dimension
        
        *Arguments*
            mesh (_mesh))
                The mesh on which the value collection is defined
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: void init(std::shared_ptr<const Mesh> mesh, std::size_t dim)
    
        Initialise MeshValueCollection with mesh and dimension
        (shared_ptr version)
        
        *Arguments*
            mesh (_mesh))
                The mesh on which the value collection is defined
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: void init(std::size_t dim)
    
        Set dimension. This function should not generally be used. It is
        for reading MeshValueCollections as the dimension is not
        generally known at construction.
        
        *Arguments*
            dim (std::size_t)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: std::size_t dim() const
    
        Return topological dimension
        
        *Returns*
            std::size_t
                The dimension.


    .. cpp:function:: std::shared_ptr<const Mesh> mesh() const
    
        Return associated mesh
        
        *Returns*
            :cpp:class:`Mesh`
                The mesh.


    .. cpp:function:: bool empty() const
    
        Return true if the subset is empty
        
        *Returns*
            bool
                True if the subset is empty.


    .. cpp:function:: std::size_t size() const
    
        Return size (number of entities in subset)
        
        *Returns*
            std::size_t
                The size.


    .. cpp:function:: bool set_value(std::size_t cell_index, std::size_t local_entity, const T& value)
    
        Set marker value for given entity defined by a cell index and
        a local entity index
        
        *Arguments*
            cell_index (std::size_t)
                The index of the cell.
            local_entity (std::size_t)
                The local index of the entity relative to the cell.
            marker_value (T)
                The value of the marker.
        
        *Returns*
            bool
                True is a new value is inserted, false if overwriting
                an existing value.


    .. cpp:function:: bool set_value(std::size_t entity_index, const T& value)
    
        Set value for given entity index
        
        *Arguments*
            entity_index (std::size_t)
                Index of the entity.
            value (T).
                The value of the marker.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            bool
                True is a new value is inserted, false if overwriting
                an existing value.


    .. cpp:function:: T get_value(std::size_t cell_index, std::size_t local_entity)
    
        Get marker value for given entity defined by a cell index and
        a local entity index
        
        *Arguments*
            cell_index (std::size_t)
                The index of the cell.
            local_entity (std::size_t)
                The local index of the entity relative to the cell.
        
        *Returns*
            marker_value (T)
                The value of the marker.


    .. cpp:function:: std::map<std::pair<std::size_t, std::size_t>, T>& values()
    
        Get all values
        
        *Returns*
            std::map<std::pair<std::size_t, std::size_t>, T>
                A map from positions to values.


    .. cpp:function:: const std::map<std::pair<std::size_t, std::size_t>, T>& values() const
    
        Get all values (const version)
        
        *Returns*
            std::map<std::pair<std::size_t, std::size_t>, T>
                A map from positions to values.


    .. cpp:function:: void clear()
    
        Clear all values


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation.


