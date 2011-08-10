
.. Documentation for the header file dolfin/quadrature/Quadrature.h

.. _programmers_reference_cpp_quadrature_quadrature:

Quadrature.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Quadrature

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: Quadrature(unsigned int n, double m=1.0)
    
        Constructor


    .. cpp:function:: int size() const
    
        Return number of quadrature points


    .. cpp:function:: double point(unsigned int i) const
    
        Return quadrature point


    .. cpp:function:: double weight(unsigned int i) const
    
        Return quadrature weight


    .. cpp:function:: double measure() const
    
        Return sum of weights (length, area, volume)


