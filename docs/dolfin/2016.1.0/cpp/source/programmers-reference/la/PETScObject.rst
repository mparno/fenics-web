
.. Documentation for the header file dolfin/la/PETScObject.h

.. _programmers_reference_cpp_la_petscobject:

PETScObject.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScObject

    This class calls SubSystemsManager to initialise PETSc.
    
    All PETSc objects must be derived from this class.


    .. cpp:function:: PETScObject()
    
        Constructor. Ensures that PETSc has been initialised.


    .. cpp:function:: static void petsc_error(int error_code, std::string filename, std::string petsc_function)
    
        Print error message for PETSc calls that return an error


