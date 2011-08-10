
.. Documentation for the header file dolfin/mesh/MeshData.h

.. _programmers_reference_cpp_mesh_meshdata:

MeshData.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshData

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    The class MeshData is a container for auxiliary mesh data,
    represented either as :cpp:class:`MeshFunction` over topological mesh
    entities, arrays or maps. Each dataset is identified by a unique
    user-specified string. Only uint-valued data are currently
    supported.
    
    Auxiliary mesh data may be attached to a mesh by users as a
    convenient way to store data associated with a mesh. It is also
    used internally by DOLFIN to communicate data associated with
    meshes. The following named mesh data are recognized by DOLFIN:
    
    Boundary indicators
    
      * "boundary_facet_cells"   -  :cpp:class:`Array` <uint> of size num_facets
      * "boundary_facet_numbers" -  :cpp:class:`Array` <uint> of size num_facets
      * "boundary_indicators"    -  :cpp:class:`Array` <uint> of size num_facets
      * "material_indicators"    -  :cpp:class:`MeshFunction` <uint> of dimension D
    
    Subdomain indicators
    
      * "cell_domains"           - :cpp:class:`MeshFunction` <uint> of dimension D
      * "interior_facet_domains" - :cpp:class:`MeshFunction` <uint> of dimension D - 1
      * "exterior_facet_domains" - :cpp:class:`MeshFunction` <uint> of dimension D - 1
    
    Facet orientation (used for assembly over interior facets)
    
      * "facet_orientation"      - :cpp:class:`MeshFunction` <uint> of dimension D - 1
    
    Sub meshes (used by the class SubMesh)
    
      * "global_vertex_indices"  - :cpp:class:`MeshFunction` <uint> of dimension 0
    
    Note to developers: use underscore in names in place of spaces.


    .. cpp:function:: MeshData(Mesh& mesh)
    
        Constructor


    .. cpp:function:: const MeshData& operator= (const MeshData& data)
    
        Assignment operator
        
        *Arguments*
            data (:cpp:class:`MeshData`)
                Another MeshData object.


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > create_mesh_function(std::string name)
    
        Create MeshFunction with given name (uninitialized)
        
        *Arguments*
            name (std::string)
                The name of the mesh function.
        
        *Returns*
            :cpp:class:`MeshFunction` <unsigned int>
                The mesh function.


    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > create_mesh_function(std::string name, uint dim)
    
        Create MeshFunction with given name and dimension
        
        *Arguments*
            name (std::string)
                The name of the mesh function.
            dim (uint)
                The dimension of the mesh function.
        
        *Returns*
            :cpp:class:`MeshFunction` <unsigned int>
                The mesh function.


    .. cpp:function:: std::vector<uint>* create_array(std::string name)
    
        Create empty array (vector) with given name
        
        *Arguments*
            name (std::string)
                The name of the array.
        
        *Returns*
            std::vector<uint>
                The array.


    .. cpp:function:: std::vector<uint>* create_array(std::string name, uint size)
    
        Create array (vector) with given name and size
        
        *Arguments*
            name (std::string)
                The name of the array.
            size (unit)
                The size (length) of the array.
        
        *Returns*
            std::vector<uint>
                The array.


    .. cpp:function:: std::map<uint, uint>* create_mapping(std::string name)
    
        Create mapping from uint to uint with given name
        
        *Arguments*
            name (std::string)
                The name of the map.
        
        *Returns*
            std::map<uint, uint>
                The map.


    .. cpp:function:: std::map<uint, std::vector<uint> >* create_vector_mapping(std::string name)
    
        Create mapping from uint to vector of uint with given name
        
        *Arguments*
            name (std::string)
                The name of the map.
        
        *Returns*
            std::map<uint, std::vector<uint> >
                The map.


    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > mesh_function(const std::string name) const
    
        Return MeshFunction with given name (returning zero if data is
        not available)
        
        *Arguments*
            name (std::string)
                The name of the MeshFunction.
        
        *Returns*
            :cpp:class:`MeshFunction` <unsigned int>
                The mesh function with given name


    .. cpp:function:: std::vector<uint>* array(const std::string name) const
    
        Return array with given name (returning zero if data is not
        available)
        
        *Arguments*
            name (std::string)
                The name of the array.
        
        *Returns*
            std::vector<uint>
                The array.


    .. cpp:function:: std::vector<uint>* array(const std::string name, uint number) const
    
        Return array with given name postfixed by " %d" (returning zero
        if data is not available)
        
        *Arguments*
            name (std::string)
                The name.
            number (uint)
                The number.
        
        *Returns*
            std::vector<uint>
                The array.


    .. cpp:function:: std::map<uint, uint>* mapping(const std::string name) const
    
        Return mapping with given name (returning zero if data is not
        available)
        
        *Arguments*
            name (std::string)
                The name of the map.
        
        *Returns*
            std::map<uint, uint>
                The map.


    .. cpp:function:: std::map<uint, std::vector<uint> >* vector_mapping(const std::string name) const
    
        Return vector mapping with given name (returning zero if data
        is not available)
        
        *Arguments*
            name (std::string)
                The name of the map
        
        *Returns*
            std::map<uint, std::vector<uint> >
                The vector mapping.


    .. cpp:function:: void erase_mesh_function(const std::string name)
    
        Erase MeshFunction with given name
        
        *Arguments*
            name (std::string)
                The name of the mesh function


    .. cpp:function:: void erase_array(const std::string name)
    
        Erase array with given name
        
        *Arguments*
            name (std::string)
                The name of the array.


    .. cpp:function:: void erase_mapping(const std::string name)
    
        Erase mapping with given name
        
        *Arguments*
            name (std::string)
                The name of the mapping.


    .. cpp:function:: void erase_vector_mapping(const std::string name)
    
        Erase vector mapping with given name
        
        *Arguments*
            name (std::string)
                The name of the vector mapping.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation.


