
.. Documentation for the header file dolfin/function/SpecialFunctions.h

.. _programmers_reference_cpp_function_specialfunctions:

SpecialFunctions.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshCoordinates

    *Parent class(es)*
    
        * :cpp:class:`Expression`
        
    This Function represents the mesh coordinates on a given mesh.


    .. cpp:function:: explicit MeshCoordinates(const Mesh& mesh)
    
        Constructor


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Evaluate function


.. cpp:class:: FacetArea

    *Parent class(es)*
    
        * :cpp:class:`Expression`
        
    This function represents the area/length of a cell facet on a given mesh.


    .. cpp:function:: explicit FacetArea(const Mesh& mesh)
    
        Constructor


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Evaluate function


