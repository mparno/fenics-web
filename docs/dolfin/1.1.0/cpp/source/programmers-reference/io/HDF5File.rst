
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
        
    .. cpp:function:: HDF5File(const std::string filename, const std::string file_mode, bool use_mpiio=true)
    
        Constructor


    .. cpp:function:: void write(const GenericVector& x, const std::string name)
    
        Write Vector to file in a format suitable for re-reading


    .. cpp:function:: void write(const Mesh& mesh, const std::string name)
    
        Write Mesh to file


    .. cpp:function:: void write(const Mesh& mesh, const std::size_t cell_dim, const std::string name)
    
        Write Mesh of given cell dimension to file in a format suitable
        for re-reading


    .. cpp:function:: void write_visualisation(const Mesh& mesh, const std::string name)
    
        Write Mesh to file for visualisation (may contain duplicate
        entities and will not preserve global indices)


    .. cpp:function:: void write_visualisation(const Mesh& mesh, const std::size_t cell_dim, const std::string name)
    
        Write Mesh of given cell dimension to file for visualisation (may
        contain duplicate entities and will not preserve global indices)


    .. cpp:function:: void read(GenericVector& x, const std::string dataset_name, const bool use_partition_from_file=true)
    
        Read vector from file


    .. cpp:function:: void read(Mesh& mesh, const std::string name)
    
        Read Mesh from file


    .. cpp:function:: bool has_dataset(const std::string dataset_name) const
    
        Check if dataset exists in HDF5 file


    .. cpp:function:: void flush()
    
        Flush buffered I/O to disk


