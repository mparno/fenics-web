
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
        
    .. cpp:function:: LocalMeshData()
    
        Create empty local mesh data


    .. cpp:function:: LocalMeshData(const Mesh& mesh)
    
        Create local mesh data for given mesh


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: void extract_mesh_data(const Mesh& mesh)
    
        Copy data from mesh


    .. cpp:function:: void broadcast_mesh_data()
    
        Broadcast mesh data from main process


    .. cpp:function:: void receive_mesh_data()
    
        Receive mesh data from main process


