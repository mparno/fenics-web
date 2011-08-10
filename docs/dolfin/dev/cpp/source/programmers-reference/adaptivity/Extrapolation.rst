
.. Documentation for the header file dolfin/adaptivity/Extrapolation.h

.. _programmers_reference_cpp_adaptivity_extrapolation:

Extrapolation.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Extrapolation

    This class implements an algorithm for extrapolating a function
    on a given function space from an approximation of that function
    on a possibly lower-order function space.
    
    This can be used to obtain a higher-order approximation of a
    computed dual solution, which is necessary when the computed
    dual approximation is in the test space of the primal problem,
    thereby being orthogonal to the residual.
    
    It is assumed that the extrapolation is computed on the same
    mesh as the original function.


    .. cpp:function:: static void extrapolate(Function& w, const Function& v)
    
        Compute extrapolation w from v


