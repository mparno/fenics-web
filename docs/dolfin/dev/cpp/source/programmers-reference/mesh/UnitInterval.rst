
.. Documentation for the header file dolfin/mesh/UnitInterval.h

.. _programmers_reference_cpp_mesh_unitinterval:

UnitInterval.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitInterval

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    A mesh of the unit interval (0, 1) with a given number of cells
    (nx) in the axial direction. The total number of intervals will
    be nx and the total number of vertices will be (nx + 1).


    .. cpp:function:: UnitInterval(uint nx=1)
    
        Create mesh of unit interval


