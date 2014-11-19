
.. Documentation for the header file dolfin/io/XDMFFile.h

.. _programmers_reference_cpp_io_xdmffile:

XDMFFile.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: XDMFFile

    *Parent class(es)*
    
        * :cpp:class:`GenericFile`
        
        * :cpp:class:`Variable`
        
    This class supports the output of meshes and functions in XDMF
    (http://www.xdmf.org) format. It creates an XML file that describes
    the data and points to a HDF5 file that stores the actual problem
    data. Output of data in parallel is supported.
    
    XDMF is not suitable for checkpointing as it may decimate
    some data.


    .. cpp:function:: explicit XDMFFile(const std::string filename)
    
        Constructor


    .. cpp:function:: void operator<< (const Mesh& mesh)
    
        Save a mesh for visualisation, with e.g. ParaView. Creates a HDF5
        file to store the mesh, and a related XDMF file with metadata.


    .. cpp:function:: void operator>> (Mesh& mesh)
    
        Read in a mesh from the associated HDF5 file


    .. cpp:function:: void operator<< (const Function& u)
    
        Save a Function to XDMF/HDF5 files for visualisation.


    .. cpp:function:: void operator<< (const std::pair<const Function*, double> ut)
    
        Save Function + time stamp to file


    .. cpp:function:: void operator<< (const MeshFunction<bool>& meshfunction)
    
        Save MeshFunction to file


    .. cpp:function:: void operator>> (MeshFunction<bool>& meshfunction)
    
        Read first MeshFunction from file


