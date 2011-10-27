
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
        


    .. cpp:function:: explicit MeshValueCollection(uint dim)
    
        Create empty mesh value collection of given dimension
        
        *Arguments*
            dim (uint)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: explicit MeshValueCollection(const MeshFunction<T>& mesh_function)
    
        Create a mesh value collection from a MeshFunction
        
        *Arguments*
            mesh_function (:cpp:class:`MeshFunction` <T>)
                The mesh function for creating a MeshValueCollection.


    .. cpp:function:: MeshValueCollection(const Mesh& mesh, const std::string filename, uint dim)
    
        Create a mesh value collection from a file.
        
        *Arguments*
            mesh (Mesh)
                A mesh associated with the collection. The mesh is used to
                map collection values to the appropriate process.
            filename (std::string)
                The XML file name.
            dim (uint)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: void set_dim(uint dim)
    
        Set the topological dimension
        
        *Arguments*
            dim (uint)
                The mesh entity dimension for the mesh value collection.


    .. cpp:function:: uint dim() const
    
        Return topological dimension
        
        *Returns*
            uint
                The dimension.


    .. cpp:function:: uint size() const
    
        Return size (number of entities in subset)
        
        *Returns*
            uint
                The size.


    .. cpp:function:: bool set_value(uint cell_index, uint local_entity, const T& value)
    
        Set marker value for given entity defined by a cell index and
        a local entity index
        
        *Arguments*
            cell_index (uint)
                The index of the cell.
            local_entity (uint)
                The local index of the entity relative to the cell.
            marker_value (T)
                The value of the marker.
        
        *Returns*
            bool
                True is a new value is inserted, false if overwriting
                an existing value.


    .. cpp:function:: bool set_value(uint entity_index, const T& value, const Mesh& mesh)
    
        Set value for given entity index
        
        *Arguments*
            entity_index (uint)
                Index of the entity.
            value (T).
                The value of the marker.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            bool
                True is a new value is inserted, false if overwriting
                an existing value.


    .. cpp:function:: std::map<std::pair<uint, uint>, T>& values()
    
        Get all values
        
        *Returns*
            std::map<std::pair<uint, uint>, T>
                A map from positions to values.


    .. cpp:function:: const std::map<std::pair<uint, uint>, T>& values() const
    
        Get all values (const version)
        
        *Returns*
            std::map<std::pair<uint, uint>, T>
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


