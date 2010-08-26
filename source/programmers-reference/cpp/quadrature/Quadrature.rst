.. Documentation for the header file dolfin/quadrature/Quadrature.h

.. _programmers_reference_cpp_quadrature_Mesh:

Quadrature.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Quadrature

    *Parent class*
    
        * :cpp:class:`Variable`
        
    .. cpp:function:: Quadrature(unsigned int n)
    
        Constructor

    .. cpp:function:: int size() const
    
        Return number of quadrature points

    .. cpp:function:: real measure() const
    
        Return sum of weights (length, area, volume)

    .. cpp:function:: real point(unsigned int i) const
    
        Return quadrature point

    .. cpp:function:: real weight(unsigned int i) const
    
        Return quadrature weight

    .. cpp:function:: virtual ~Quadrature()
    
        Destructor

