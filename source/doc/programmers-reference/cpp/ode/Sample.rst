.. Documentation for the header file dolfin/ode/Sample.h

.. _programmers_reference_cpp_ode_sample:

Sample.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Sample

    *Parent class*
    
        * :cpp:class:`Variable`
        
    Sample of solution values at a given point.

    .. cpp:function:: Sample(TimeSlab& timeslab, real t, std::string name, std::string label)
    
        Constructor

    .. cpp:function:: uint size() const
    
        Return number of components

    .. cpp:function:: real t() const
    
        Return time t

    .. cpp:function:: real u(uint index) const
    
        Return value of component with given index

    .. cpp:function:: real k(uint index) const
    
        Return time step for component with given index

    .. cpp:function:: real r(uint index) const
    
        Return residual for component with given index

