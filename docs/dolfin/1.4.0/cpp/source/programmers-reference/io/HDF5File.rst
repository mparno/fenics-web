
.. Documentation for the header file dolfin/io/HDF5File.h

.. _programmers_reference_cpp_io_hdf5file:

HDF5File.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: HDF5File

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: HDF5File(MPI_Comm comm, const std::string filename, const std::string file_mode)
    
        Constructor. file_mode should "a" (append), "w" (write) or "r"
        (read).


    .. cpp:function:: void close()
    
        Close file


    .. cpp:function:: void write(const std::vector<Point>& points, const std::string name)
    
        Write points to file


    .. cpp:function:: void write(const std::vector<double>& values, const std::string name)
    
        Write simple vector of double to file


    .. cpp:function:: void write(const GenericVector& x, const std::string name)
    
        Write Vector to file in a format suitable for re-reading


    .. cpp:function:: void read(GenericVector& x, const std::string dataset_name, const bool use_partition_from_file) const
    
        Read vector from file and optionally re-use any partitioning
        that is available in the file


    .. cpp:function:: void write(const Mesh& mesh, const std::string name)
    
        Write Mesh to file in a format suitable for re-reading


    .. cpp:function:: void write(const Mesh& mesh, const std::size_t cell_dim, const std::string name)
    
        Write Mesh of given cell dimension to file in a format
        suitable for re-reading


    .. cpp:function:: void write(const Function& u, const std::string name)
    
        Write Function to file in a format suitable for re-reading


    .. cpp:function:: void write(const Function& u, const std::string name, double timestamp)
    
        Write Function to file with a timestamp


    .. cpp:function:: void read(Function& u, const std::string name)
    
        Read Function from file and distribute data according to
        the Mesh and dofmap associated with the Function


    .. cpp:function:: void read(Mesh& mesh, const std::string name, bool use_partition_from_file) const
    
        Read Mesh from file and optionally re-use any partition data
        in the file


    .. cpp:function:: void write(const MeshFunction<std::size_t>& meshfunction, const std::string name)
    
        Write MeshFunction to file in a format suitable for re-reading


    .. cpp:function:: void write(const MeshFunction<int>& meshfunction, const std::string name)
    
        Write MeshFunction to file in a format suitable for re-reading


    .. cpp:function:: void write(const MeshFunction<double>& meshfunction, const std::string name)
    
        Write MeshFunction to file in a format suitable for re-reading


    .. cpp:function:: void write(const MeshFunction<bool>& meshfunction, const std::string name)
    
        Write MeshFunction to file in a format suitable for re-reading


    .. cpp:function:: void read(MeshFunction<std::size_t>& meshfunction, const std::string name) const
    
        Read MeshFunction from file


    .. cpp:function:: void read(MeshFunction<int>& meshfunction, const std::string name) const
    
        Read MeshFunction from file


    .. cpp:function:: void read(MeshFunction<double>& meshfunction, const std::string name) const
    
        Read MeshFunction from file


    .. cpp:function:: void read(MeshFunction<bool>& meshfunction, const std::string name) const
    
        Read MeshFunction from file


    .. cpp:function:: void write(const MeshValueCollection<std::size_t>& mesh_values, const std::string name)
    
        Write MeshValueCollection to file


    .. cpp:function:: void write(const MeshValueCollection<double>& mesh_values, const std::string name)
    
        Write MeshValueCollection to file


    .. cpp:function:: void write(const MeshValueCollection<bool>& mesh_values, const std::string name)
    
        Write MeshValueCollection to file


    .. cpp:function:: void read(MeshValueCollection<std::size_t>& mesh_values, const std::string name) const
    
        Read MeshValueCollection from file


    .. cpp:function:: void read(MeshValueCollection<double>& mesh_values, const std::string name) const
    
        Read MeshValueCollection from file


    .. cpp:function:: void read(MeshValueCollection<bool>& mesh_values, const std::string name) const
    
        Read MeshValueCollection from file


    .. cpp:function:: bool has_dataset(const std::string dataset_name) const
    
        Check if dataset exists in HDF5 file


    .. cpp:function:: void flush()
    
        Flush buffered I/O to disk


