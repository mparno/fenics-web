
.. Documentation for the header file dolfin/mesh/LocalMeshData.h

.. _programmers_reference_cpp_mesh_localmeshdata:

LocalMeshData.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalMeshData

    *Parent class(es)*
    
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


    .. cpp:function:: LocalMeshData()
    
        Create empty local mesh data


    .. cpp:function:: LocalMeshData(const Mesh& mesh)
    
        Create local mesh data for given mesh


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


