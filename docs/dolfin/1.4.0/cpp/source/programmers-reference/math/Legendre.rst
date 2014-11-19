
.. Documentation for the header file dolfin/math/Legendre.h

.. _programmers_reference_cpp_math_legendre:

Legendre.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Legendre

    Interface for computing Legendre polynomials via Boost.


    .. cpp:function:: static double eval(std::size_t n, double x)
    
        Evaluate polynomial of order n at point x


    .. cpp:function:: static double ddx(std::size_t n, double x)
    
        Evaluate first derivative of polynomial of order n at point x


    .. cpp:function:: static double d2dx(std::size_t n, double x)
    
        Evaluate second derivative of polynomial of order n at point x


