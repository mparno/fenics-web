
.. Documentation for the header file dolfin/adaptivity/TimeSeries.h

.. _programmers_reference_cpp_adaptivity_timeseries:

TimeSeries.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TimeSeries

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class stores a time series of objects to file(s) in a
    binary format which is efficient for reading and writing.
    
    When objects are retrieved, the object stored at the time
    closest to the given time will be used.
    
    A new time series will check if values have been stored to
    file before (for a series with the same name) and in that
    case reuse those values. If new values are stored, old
    values will be cleared.


    .. cpp:function:: TimeSeries(std::string name)
    
        Create empty time series
        
        *Arguments*
            name (std::string)
                The time series name


    .. cpp:function:: void store(const GenericVector& vector, double t)
    
        Store vector at given time
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector to be stored.
            t (double)
                The time.


    .. cpp:function:: void store(const Mesh& mesh, double t)
    
        Store mesh at given time
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to be stored.
            t (double)
                The time.


    .. cpp:function:: void retrieve(GenericVector& vector, double t, bool interpolate=true) const
    
        Retrieve vector at given time
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector (values to be retrieved).
            t (double)
                The time.
            interpolate (bool)
                Optional argument: If true (default), interpolate
                time samples closest to t if t is not present.


    .. cpp:function:: void retrieve(Mesh& mesh, double t) const
    
        Retrieve mesh at given time
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh (values to be retrieved).
            t (double)
                The time.


    .. cpp:function:: Array<double> vector_times() const
    
        Return array of sample times for vectors
        
        *Returns*
            :cpp:class:`Array` <double>
                The times.


    .. cpp:function:: Array<double> mesh_times() const
    
        Return array of sample times for meshes
        
        *Returns*
            :cpp:class:`Array` <double>
                The times.


    .. cpp:function:: void clear()
    
        Clear time series


    .. cpp:function:: static std::string filename_data(std::string series_name, std::string type_name, uint index)
    
        Return filename for data
        
        *Arguments*
            series_name (std::string)
                The time series name
            type_name (std::string)
                The type of data
            index (uint)
                The index
        
        *Returns*
            std::string
                The filename


    .. cpp:function:: static std::string filename_times(std::string series_name, std::string type_name)
    
        Return filename for times
        
        *Arguments*
            series_name (std::string)
                The time series name
            type_name (std::string)
                The type of data
        
        *Returns*
            std::string
                The filename


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


