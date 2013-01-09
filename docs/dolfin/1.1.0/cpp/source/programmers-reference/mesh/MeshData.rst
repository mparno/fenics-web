
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
    user-specified string. Only std::size_t-valued data are currently
    supported.
    
    Auxiliary mesh data may be attached to a mesh by users as a
    convenient way to store data associated with a mesh. It is also
    used internally by DOLFIN to communicate data associated with
    meshes. The following named mesh data are recognized by DOLFIN:
    
    Facet orientation (used for assembly over interior facets)
    
      * "facet_orientation"     - :cpp:class:`MeshFunction` <std::size_t> of dimension D - 1
    
    Sub meshes (used by the class SubMesh)
    
      * "parent_vertex_indices" - :cpp:class:`MeshFunction` <std::size_t> of dimension 0
    
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


    .. cpp:function:: boost::shared_ptr<MeshFunction<std::size_t> > create_mesh_function(std::string name)
    
        Create MeshFunction with given name (uninitialized)
        
        *Arguments*
            name (std::string)
                The name of the mesh function.
        
        *Returns*
            :cpp:class:`MeshFunction` <std::size_t>
                The mesh function.


    .. cpp:function:: boost::shared_ptr<MeshFunction<std::size_t> > create_mesh_function(std::string name, std::size_t dim)
    
        Create MeshFunction with given name and dimension
        
        *Arguments*
            name (std::string)
                The name of the mesh function.
            dim (std::size_t)
                The dimension of the mesh function.
        
        *Returns*
            :cpp:class:`MeshFunction` <std::size_t>
                The mesh function.


    .. cpp:function:: boost::shared_ptr<std::vector<std::size_t> > create_array(std::string name)
    
        Create empty array (vector) with given name
        
        *Arguments*
            name (std::string)
                The name of the array.
        
        *Returns*
            std::vector<std::size_t>
                The array.


    .. cpp:function:: boost::shared_ptr<std::vector<std::size_t> > create_array(std::string name, std::size_t size)
    
        Create array (vector) with given name and size
        
        *Arguments*
            name (std::string)
                The name of the array.
            size (std::size_t)
                The size (length) of the array.
        
        *Returns*
            std::vector<std::size_t>
                The array.


    .. cpp:function:: boost::shared_ptr<MeshFunction<std::size_t> > mesh_function(const std::string name) const
    
        Return MeshFunction with given name (returning zero if data is
        not available)
        
        *Arguments*
            name (std::string)
                The name of the MeshFunction.
        
        *Returns*
            :cpp:class:`MeshFunction` <std::size_t>
                The mesh function with given name


    .. cpp:function:: boost::shared_ptr<std::vector<std::size_t> > array(const std::string name) const
    
        Return array with given name (returning zero if data is not
        available)
        
        *Arguments*
            name (std::string)
                The name of the array.
        
        *Returns*
            std::vector<std::size_t>
                The array.


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


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation.


