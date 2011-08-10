
.. Documentation for the header file dolfin/parameter/GlobalParameters.h

.. _programmers_reference_cpp_parameter_globalparameters:

GlobalParameters.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GlobalParameters

    *Parent class(es)*
    
        * :cpp:class:`Parameters`
        
    This class defines the global DOLFIN parameter database.


    .. cpp:function:: GlobalParameters()
    
        Constructor


    .. cpp:function:: void parse(int argc, char* argv[])
    
        Parse parameters from command-line


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


