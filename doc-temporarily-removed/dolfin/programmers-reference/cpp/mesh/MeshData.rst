.. Documentation for the header file dolfin/mesh/MeshData.h

.. _programmers_reference_cpp_mesh_meshdata:

MeshData.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MeshData

    *Parent class*
    
        * :cpp:class:`Variable`
        
    The class MeshData is a container for auxiliary mesh data,
    represented either as MeshFunctions over topological mesh
    entities, arrays or maps. Each dataset is identified by a unique
    user-specified string. Only uint-valued data are currently
    supported.
    
    The following named mesh data are recognized by DOLFIN:
    
    Boundary indicators
    
      "boundary_facet_cells"   - Array<uint> of size num_facets
      "boundary_facet_numbers" - Array<uint> of size num_facets
      "boundary_indicators"    - Array<uint> of size num_facets
      "material_indicators"    - MeshFunction<uint> of dimension D
    
    Subdomain indicators
    
      "cell_domains"           - MeshFunction<uint> of dimension D
      "interior_facet_domains" - MeshFunction<uint> of dimension D - 1
      "exterior_facet_domains" - MeshFunction<uint> of dimension D - 1
    
    Facet orientation (used for assembly over interior facets)
    
      "facet orientation" - MeshFunction<uint> of dimension D - 1
    
    Boundary extraction
    
      "vertex map" - MeshFunction<uint> of dimension 0
      "cell map"   - MeshFunction<uint> of dimension D
    
    Mesh partitioning
    
      (moved to ParallelData) "global entity indices %d" - MeshFunction<uint> of dimension 0, 1, ..., D
      (moved to ParallelData) "exterior facets"          - MeshFunction<uint> of dimension D - 1
      (moved to ParallelData) "num global entities"      - Array<uint> of size D + 1
      (moved to ParallelData) "overlap"                  - vector mapping
    
    Sub meshes
    
      "global vertex indices" - MeshFunction<uint> of dimension 0
    
    Mesh coloring
    
      "colors-%D-%d-%1"   - MeshFunction<uint> of dimension D with colors based on connectivity %d
      "num colored cells" - Array<uint> listing the number of cells of each color
      "colored cells %d"  - Array<uint> of cell indices with colors 0, 1, 2, ...

    .. cpp:function:: MeshData(Mesh& mesh)
    
        Constructor

    .. cpp:function:: const MeshData& operator= (const MeshData& data)
    
        Assignment operator

    .. cpp:function:: void clear()
    
        Clear all data

    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > create_mesh_function(std::string name)
    
        Create MeshFunction with given name (uninitialized)

    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > create_mesh_function(std::string name, uint dim)
    
        Create MeshFunction with given name and dimension

    .. cpp:function:: std::vector<uint>* create_array(std::string name)
    
        Create empty array (vector) with given name

    .. cpp:function:: std::vector<uint>* create_array(std::string name, uint size)
    
        Create array (vector) with given name and size

    .. cpp:function:: std::map<uint, uint>* create_mapping(std::string name)
    
        Create mapping from uint to uint with given name

    .. cpp:function:: std::map<uint, std::vector<uint> >* create_vector_mapping(std::string name)
    
        Create mapping from uint to vector of uint with given name

    .. cpp:function:: boost::shared_ptr<MeshFunction<unsigned int> > mesh_function(const std::string name) const
    
        Return MeshFunction with given name (returning zero if data is not available)

    .. cpp:function:: std::vector<uint>* array(const std::string name) const
    
        Return array with given name (returning zero if data is not available)

    .. cpp:function:: std::vector<uint>* array(const std::string name, uint number) const
    
        Return array with given name postfixed by " %d" (returning zero if data is not available)

    .. cpp:function:: std::map<uint, uint>* mapping(const std::string name) const
    
        Return mapping with given name (returning zero if data is not available)

    .. cpp:function:: std::map<uint, std::vector<uint> >* vector_mapping(const std::string name) const
    
        Return vector mapping with given name (returning zero if data is not available)

    .. cpp:function:: void erase_mesh_function(const std::string name)
    
        Erase MeshFunction with given name

    .. cpp:function:: void erase_array(const std::string name)
    
        Erase array with given name

    .. cpp:function:: void erase_mapping(const std::string name)
    
        Erase mapping with given name

    .. cpp:function:: void erase_vector_mapping(const std::string name)
    
        Erase vector mapping with given name

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

