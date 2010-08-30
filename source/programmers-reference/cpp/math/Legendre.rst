.. Documentation for the header file dolfin/math/Legendre.h

.. _programmers_reference_cpp_math_legendre:

Legendre.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Legendre

    Legendre polynomial of given degree n on the interval [-1,1].
    
      P0(x) = 1
      P1(x) = x
      P2(x) = (3x^2 - 1) / 2
      ...
    
    The function values and derivatives are computed using
    three-term recurrence formulas.

    .. cpp:function:: real d2dx(real x)
    
        Evaluation of second derivative at given point

    .. cpp:function:: real ddx(real x)
    
        Evaluation of derivative at given point

    .. cpp:function:: real eval(uint nn, real x)
    
        Evaluation of arbitrary order, nn <= n (useful ie in RadauQuadrature)

    .. cpp:function:: real operator() (real x)
    
        Evaluation at given point

