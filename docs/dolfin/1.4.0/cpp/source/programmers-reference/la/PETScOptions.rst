
.. Documentation for the header file dolfin/la/PETScOptions.h

.. _programmers_reference_cpp_la_petscoptions:

PETScOptions.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScOptions

    These class provides static functions that permit users to set
    and retreive PETSc options via the PETSc option/parameter
    system. The option should be prefixed by '-', e.g.
    
        PETScOptions::set("mat_mumps_icntl_14", 40.0);


    .. cpp:function:: static void set(std::string option)
    
        Set PETSc option that takes no value


    .. cpp:function:: static void set(std::string option, bool value)
    
        Set PETSc boolean option


    .. cpp:function:: static void set(std::string option, int value)
    
        Set PETSc integer option


    .. cpp:function:: static void set(std::string option, double value)
    
        Set PETSc double option


    .. cpp:function:: static void set(std::string option, std::string value)
    
        Set PETSc string option


    .. cpp:function:: static void set(std::string option, const T value)
    
        Genetic function for setting PETSc option


    .. cpp:function:: static void clear(std::string option)
    
        Clear PETSc option


