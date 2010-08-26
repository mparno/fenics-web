.. Documentation for the header file dolfin/mesh/LocalMeshData.h

.. _programmers_reference_cpp_mesh_localmeshdata:

LocalMeshData.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: LocalMeshData

    *Parent class*
    
        * :cpp:class:`Variable`
        
        This class stores mesh data on a local processor corresponding
        to a portion of a (larger) global mesh.
        
        Note that the data stored in this class does typically not
        correspond to a topologically connected mesh; it merely stores a
        list of vertex coordinates, a list of cell-vertex mappings and a
        list of global vertex numbers for the locally stored vertices.
        
        It is typically used for parsing meshes in parallel from mesh
        XML files. After local mesh data has been parsed on each
        processor, a subsequent repartitioning takes place: first a
        geometric partitioning of the vertices followed by a
        redistribution of vertex and cell data, and then a topological
        partitioning again followed by redistribution of vertex and cell
        data, at that point corresponding to topologically connected
        meshes instead of local mesh data.

    .. cpp:function:: // FIXME: Provide a better public interface rather than using 'friend class'
                       
                       LocalMeshData()
    
        Create empty local mesh data

    .. cpp:function:: LocalMeshData(const Mesh& mesh)
    
        Create local mesh data for given mesh

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: std::vector<std::vector<double> > vertex_coordinates
    
        Coordinates for all vertices stored on local processor

    .. cpp:function:: std::vector<std::vector<uint> > cell_vertices
    
        Global vertex indices for all cells stored on local processor

    .. cpp:function:: std::vector<uint> global_cell_indices
    
        Global cell numbers for all cells stored on local processor

    .. cpp:function:: std::vector<uint> vertex_indices
    
        Global vertex indices for all vertices stored on local processor

    .. cpp:function:: typedef XMLLocalMeshData XMLHandler
    
        Define XMLHandler for use in new XML reader/writer

    .. cpp:function:: uint gdim
    
        Geometrical dimension

    .. cpp:function:: uint num_global_cells
    
        Global number of cells

    .. cpp:function:: uint num_global_vertices
    
        Global number of vertices

    .. cpp:function:: uint tdim
    
        Topological dimension

    .. cpp:function:: void broadcast_mesh_data()
    
        Broadcast mesh data from main process

    .. cpp:function:: void clear()
    
        Clear all data

    .. cpp:function:: void extract_mesh_data(const Mesh& mesh)
    
        Copy data from mesh

    .. cpp:function:: void receive_mesh_data()
    
        Receive mesh data from main process

    .. cpp:function:: ~LocalMeshData()
    
        Destructor

