.. Documentation for the header file dolfin/quadrature/Quadrature.h

.. _programmers_reference_cpp_quadrature_quadrature:

Quadrature.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Quadrature

    *Parent class*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: Quadrature(unsigned int n, real m=1.0)
    
        Constructor

    .. cpp:function:: int size() const
    
        Return number of quadrature points

    .. cpp:function:: real point(unsigned int i) const
    
        Return quadrature point

    .. cpp:function:: real weight(unsigned int i) const
    
        Return quadrature weight

    .. cpp:function:: real measure() const
    
        Return sum of weights (length, area, volume)

