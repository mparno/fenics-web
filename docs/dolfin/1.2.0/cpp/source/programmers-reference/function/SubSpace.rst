
.. Documentation for the header file dolfin/function/SubSpace.h

.. _programmers_reference_cpp_function_subspace:

SubSpace.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SubSpace

    *Parent class(es)*
    
        * :cpp:class:`FunctionSpace`
        
    This class represents a subspace (component) of a function space.
    
    The subspace is specified by an array of indices. For example,
    the array [3, 0, 2] specifies subspace 2 of subspace 0 of
    subspace 3.
    
    A typical example is the function space W = V x P for Stokes.
    Here, V = W[0] is the subspace for the velocity component and
    P = W[1] is the subspace for the pressure component. Furthermore,
    W[0][0] = V[0] is the first component of the velocity space etc.


    .. cpp:function:: SubSpace(const FunctionSpace& V, std::size_t component)
    
        Create subspace for given component (one level)


    .. cpp:function:: SubSpace(const FunctionSpace& V, std::size_t component, std::size_t sub_component)
    
        Create subspace for given component (two levels)


    .. cpp:function:: SubSpace(const FunctionSpace& V, const std::vector<std::size_t>& component)
    
        Create subspace for given component (n levels)


