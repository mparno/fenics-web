.. Documentation for the header file dolfin/function/SpecialFunctions.h

.. _programmers_reference_cpp_function_specialfunctions:

SpecialFunctions.h
==================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: MeshCoordinates

    *Parent class*
    
        * :cpp:class:`Expression`
        
        This Function represents the mesh coordinates on a given mesh.

    .. cpp:function:: MeshCoordinates(const Mesh& mesh)
    
        Constructor

    .. cpp:function:: void eval(Array<double>& values, const Data& data) const
    
        Evaluate function

.. cpp:class:: CellSize

    *Parent class*
    
        * :cpp:class:`Expression`
        
        This Function represents the local cell size on a given mesh.

    .. cpp:function:: CellSize(const Mesh& mesh)
    
        Constructor

    .. cpp:function:: void eval(Array<double>& values, const Data& data) const
    
        Evaluate function

.. cpp:class:: FacetArea

    *Parent class*
    
        * :cpp:class:`Expression`
        
        This function represents the area/length of a cell facet on a given mesh.

    .. cpp:function:: FacetArea(const Mesh& mesh)
    
        Constructor

    .. cpp:function:: void eval(Array<double>& values, const Data& data) const
    
        Evaluate function

