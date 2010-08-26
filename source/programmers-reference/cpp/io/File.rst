.. Documentation for the header file dolfin/io/File.h

.. _programmers_reference_cpp_io_Mesh:

File.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: File

        A File represents a data file for reading and writing objects.
        Unless specified explicitly, the format is determined by the
        file name suffix.
        A list of objects that can be read/written to file can be found in
        GenericFile.h

    .. cpp:function:: File(const std::string filename, Type type, std::string encoding = "ascii")
    
        Create a file with given name and type (format)

    .. cpp:function:: File(const std::string filename, std::string encoding = "ascii")
    
        Create a file with given name

    .. cpp:function:: File(std::ostream& outstream)
    
        Create a outfile object writing to stream

    .. cpp:function:: enum Type
    
        File formats

    .. cpp:function:: static bool exists(std::string filename)
    
        Check if file exists

    .. cpp:function:: template<class T> void operator<<(const T& t)
    
        Write to file

    .. cpp:function:: template<class T> void operator>>(T& t)
    
        Read from file

    .. cpp:function:: void operator<<(const Function& u)
    
        Write Function to file

    .. cpp:function:: void operator<<(const std::pair<const Function*, double> u)
    
        Write Function to file (with time)

    .. cpp:function:: ~File()
    
        Destructor

