
.. Documentation for the header file dolfin/plot/FunctionPlotData.h

.. _programmers_reference_cpp_plot_functionplotdata:

FunctionPlotData.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FunctionPlotData

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class is used for communicating plot data for functions
    to and from (XML) files. It is used by DOLFIN for plotting
    Function objects. The data is stored as a mesh and a vector
    of interpolated vertex values.


    .. cpp:function:: FunctionPlotData(const GenericFunction& v, const Mesh& mesh)
    
        Create plot data for given function


    .. cpp:function:: FunctionPlotData()
    
        Create empty data to be read from file


    .. cpp:function:: GenericVector& vertex_values() const
    
        Return vertex values


