.. Documentation for the header file dolfin/function/Data.h

.. _programmers_reference_cpp_function_Mesh:

Data.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Data

        This class holds data for function evaluation, including
        the coordinates x, the time t, and auxiliary data that a
        function may depend on.

    .. cpp:function:: Data()
    
        Constructor

    .. cpp:function:: Point normal() const
    
        Return current facet normal (if available)

    .. cpp:function:: bool on_facet() const
    
        Check if we are on a facet

    .. cpp:function:: const Array<double> x
    
        The coordinates

    .. cpp:function:: const Cell& cell() const
    
        Return current cell (if available)

    .. cpp:function:: const ufc::cell& ufc_cell() const
    
        Return current UFC cell (if available)

    .. cpp:function:: uint facet() const
    
        Return current facet (if available)

    .. cpp:function:: uint geometric_dimension() const
    
        Return geometric dimension of cell

    .. cpp:function:: void clear()
    
        Clear all cell data

    .. cpp:function:: void set(const Cell& dolfin_cell, const ufc::cell& ufc_cell, int local_facet)
    
        Set cell and facet data

    .. cpp:function:: void set(const ufc::cell& ufc_cell, const double* x)
    
        Set UFC cell and coordinate

    .. cpp:function:: ~Data()
    
        Destructor

