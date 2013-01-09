
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


    .. cpp:function:: UnitIntervalMesh(std::size_t nx=1)
    
        Create mesh of unit interval


