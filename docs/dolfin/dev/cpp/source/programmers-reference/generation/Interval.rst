
.. Documentation for the header file dolfin/generation/Interval.h

.. _programmers_reference_cpp_generation_interval:

Interval.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Interval

    *Parent class(es)*
    
        * :cpp:class:`IntervalMesh`
        
    Interval mesh of the 1D line [a,b].  Given the number of cells
    (nx) in the axial direction, the total number of intervals will
    be nx and the total number of vertices will be (nx + 1).
    
    This class is deprecated. Use :cpp:class:`IntervalMesh`.


    .. cpp:function:: Interval(std::size_t nx, double a, double b)
    
        Constructor
        
        *Arguments*
            nx (std::size_t)
                The number of cells.
            a (double)
                The minimum point (inclusive).
            b (double)
                The maximum point (inclusive).
        
        *Example*
            .. code-block:: c++
        
                // Create a mesh of 25 cells in the interval [-1,1]
                Interval mesh(25, -1.0, 1.0);
        


