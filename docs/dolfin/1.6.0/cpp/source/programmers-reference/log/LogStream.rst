
.. Documentation for the header file dolfin/log/LogStream.h

.. _programmers_reference_cpp_log_logstream:

LogStream.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LogStream

    This class provides functionality similar to standard C++
    streams (std::cout, std::endl) for output but working through
    the DOLFIN log system.


    .. cpp:function:: enum Type
    
        Stream types


    .. cpp:function:: LogStream(Type type)
    
        Create log stream of given type


    .. cpp:function:: LogStream& operator<< (const LogStream& stream)
    
        Output for log stream


    .. cpp:function:: LogStream& operator<< (const std::string& s)
    
        Output for string


    .. cpp:function:: LogStream& operator<< (int a)
    
        Output for int


    .. cpp:function:: LogStream& operator<< (unsigned int a)
    
        Output for unsigned int


    .. cpp:function:: LogStream& operator<< (long int a)
    
        Output for long int


    .. cpp:function:: LogStream& operator<< (long unsigned int a)
    
        Output for long int


    .. cpp:function:: LogStream& operator<< (double a)
    
        Output for double


    .. cpp:function:: LogStream& operator<< (std::complex<double> z)
    
        Output for std::complex<double>


    .. cpp:function:: LogStream& operator<< (const Variable& variable)
    
        Output for variable (calling str() method)


    .. cpp:function:: LogStream& operator<< (const MeshEntity& entity)
    
        Output for mesh entity (not subclass of Variable for efficiency)


    .. cpp:function:: LogStream& operator<< (const Point& point)
    
        Output for point (not subclass of Variable for efficiency)


