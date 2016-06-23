
.. Documentation for the header file dolfin/function/LagrangeInterpolator.h

.. _programmers_reference_cpp_function_lagrangeinterpolator:

LagrangeInterpolator.h
======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LagrangeInterpolator

    This class interpolates efficiently from a GenericFunction to a
    Lagrange Function


    .. cpp:function:: static void interpolate(Function& u, const Expression& u0)
    
        Interpolate Expression
        
        *Arguments*
            u  (:cpp:class:`Function`)
                The resulting Function
            u0 (:cpp:class:`Expression`)
                The Expression to be interpolated.


    .. cpp:function:: static void interpolate(Function& u, const Function& u0)
    
        Interpolate function (on possibly non-matching meshes)
        
        *Arguments*
            u  (:cpp:class:`Function`)
                The resulting Function
            u0 (:cpp:class:`Function`)
                The Function to be interpolated.


