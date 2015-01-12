
.. Documentation for the header file dolfin/generation/UnitIntervalMesh.h

.. _programmers_reference_cpp_generation_unitintervalmesh:

UnitIntervalMesh.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitIntervalMesh

    *Parent class(es)*
    
        * :cpp:class:`IntervalMesh`
        
    A mesh of the unit interval (0, 1) with a given number of cells
    (nx) in the axial direction. The total number of intervals will
    be nx and the total number of vertices will be (nx + 1).


    .. cpp:function:: UnitIntervalMesh(std::size_t nx)
    
        Constructor
        
        *Arguments*
            nx (std::size_t)
                The number of cells.
        
        *Example*
            .. code-block:: c++
        
                // Create a mesh of 25 cells in the interval [0,1]
                UnitIntervalMesh mesh(25);
        


    .. cpp:function:: UnitIntervalMesh(MPI_Comm comm, std::size_t nx)
    
        Constructor
        
        *Arguments*
            comm (MPI_Comm)
                MPI communicator
            nx (std::size_t)
                The number of cells.
        
        *Example*
            .. code-block:: c++
        
                // Create a mesh of 25 cells in the interval [0,1]
                UnitIntervalMesh mesh(MPI_COMM_WORLD, 25);
        


