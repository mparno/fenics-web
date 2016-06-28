
.. Documentation for the header file dolfin/io/XDMFFile.h

.. _programmers_reference_cpp_io_xdmffile:

XDMFFile.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: XDMFFile

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class supports the output of meshes and functions in XDMF
    (http://www.xdmf.org) format. It creates an XML file that
    describes the data and points to a HDF5 file that stores the
    actual problem data. Output of data in parallel is supported.
    
    XDMF is not suitable for checkpointing as it may decimate some
    data.


    .. cpp:function:: enum class Encoding
    
        File encoding type


    .. cpp:function:: enum class UseFilePartition
    
        Re-use any partition stored in file


    .. cpp:function:: XDMFFile(const std::string filename)
    
        Constructor


    .. cpp:function:: XDMFFile(MPI_Comm comm, const std::string filename)
    
        Constructor


    .. cpp:function:: void write(const Mesh& mesh, Encoding encoding=Encoding::HDF5)
    
        Save a mesh to XDMF format, either using an associated HDF5
        file, or storing the data inline as XML Create function on
        given function space
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                A mesh to save.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const Function& u, Encoding encoding=Encoding::HDF5)
    
        Save a Function to XDMF file for visualisation, using an
        associated HDF5 file, or storing the data inline as XML.
        
        *Arguments*
            u (:cpp:class:`Function`)
                A function to save.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const Function& u, double t, Encoding encoding=Encoding::HDF5)
    
        Save a Function with timestamp to XDMF file for visualisation,
        using an associated HDF5 file, or storing the data inline as
        XML.
        
        *Arguments*
            u (:cpp:class:`Function`)
                A function to save.
            t (_double_)
                Timestep
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const MeshFunction<bool>& meshfunction, Encoding encoding=Encoding::HDF5)
    
        Save MeshFunction to file using an associated HDF5 file, or
        storing the data inline as XML.
        
        *Arguments*
            meshfunction (:cpp:class:`MeshFunction`)
                A meshfunction to save.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const MeshValueCollection<std::size_t>& mvc, Encoding encoding=Encoding::HDF5)
    
        Write out mesh value collection (subset) using an associated
        HDF5 file, or storing the data inline as XML.
        
        *Arguments*
            mvc (_MeshValueCollection<std::size_t>_)
                A list of points to save.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const std::vector<Point>& points, Encoding encoding=Encoding::HDF5)
    
        Save a cloud of points to file using an associated HDF5 file,
        or storing the data inline as XML.
        
        *Arguments*
            points (_std::vector<Point>_)
                A list of points to save.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void write(const std::vector<Point>& points, const std::vector<double>& values, Encoding encoding=Encoding::HDF5)
    
        Save a cloud of points, with scalar values using an associated
        HDF5 file, or storing the data inline as XML.
        
        *Arguments*
            points (_std::vector<Point>_)
                A list of points to save.
            values (_std::vector<double>_)
                A list of values at each point.
            encoding (_Encoding_)
                Encoding to use: HDF5 or ASCII
        


    .. cpp:function:: void read(Mesh& mesh, UseFilePartition use_file_partition=UseFilePartition::no)
    
        Read in a mesh from the associated HDF5 file, optionally using
        stored partitioning, if possible when the same number of
        processes are being used.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
        
            use_partition_from_file (_UseFilePartition_)
                Use the existing partition information in HDF5 file
        


    .. cpp:function:: void read(MeshFunction<bool>& meshfunction)
    
        Read first MeshFunction from file


