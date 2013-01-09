
.. Documentation for the header file dolfin/io/File.h

.. _programmers_reference_cpp_io_file:

File.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: File

    A File represents a data file for reading and writing objects.
    Unless specified explicitly, the format is determined by the
    file name suffix.
    A list of objects that can be read/written to file can be found in
    GenericFile.h. Compatible file formats include:
        * XDMF   (.xdmf)
        * XML    (.xml)
        * VTK    (.pvd)
        * RAW    (.raw)
        * XYZ    (.xyz)
        * Binary (.bin)
        * SVG    (.svg)


    .. cpp:function:: File(const std::string filename, std::string encoding="ascii")
    
        Create a file with given name
        
        *Arguments*
            filename (std::string)
                Name of file.
            encoding (std::string)
                Optional argument specifying encoding, ascii is default.
        
        *Example*
           .. code-block:: c++
        
                // Save solution to file
                File file("solution.pvd");
                file << u;
        
                // Read mesh data from file
                File mesh_file("mesh.xml");
                mesh_file >> mesh;
        
                // Using compressed binary format
                File comp_file("solution.pvd", "compressed");
        


    .. cpp:function:: File(const std::string filename, Type type, std::string encoding="ascii")
    
        Create a file with given name and type (format)
        
        *Arguments*
            filename (std::string)
                Name of file.
            type (Type)
                File format.
            encoding (std::string)
                Optional argument specifying encoding, ascii is default.
        
        *Example*
            .. code-block:: c++
        
                File file("solution", vtk);
        


    .. cpp:function:: File(std::ostream& outstream)
    
        Create an outfile object writing to stream
        
        *Arguments*
            outstream (std::ostream)
                The stream.


    .. cpp:function:: void operator>>(T& t)
    
        Read from file


    .. cpp:function:: void operator<<(const std::pair<const Function*, double> u)
    
        Write Function to file with time
        
        *Example*
            .. code-block:: c++
        
                File file("solution.pvd", "compressed");
                file << std::make_pair<const Function*, double>(&u, t);
        


    .. cpp:function:: void operator<<(const T& t)
    
        Write object to file


    .. cpp:function:: static bool exists(std::string filename)
    
        Check if file exists
        
        *Arguments*
            filename (std::string)
                Name of file.
        
        *Returns*
            bool
                True if the file exists.


    .. cpp:function:: static void create_parent_path(std::string filename)
    
        
        *Arguments*
            filename (std::string)
                Name of file / path.


