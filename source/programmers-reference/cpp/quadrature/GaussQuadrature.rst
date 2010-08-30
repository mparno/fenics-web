.. Documentation for the header file dolfin/quadrature/GaussQuadrature.h

.. _programmers_reference_cpp_quadrature_gaussquadrature:

GaussQuadrature.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GaussQuadrature

    *Parent class*
    
        * :cpp:class:`GaussianQuadrature`
        
    Gauss (Gauss-Legendre) quadrature on the interval [-1,1].
    The n quadrature points are given by the zeros of the
    n:th Legendre Pn(x).
    
    The quadrature points are computed using Newton's method, and
    the quadrature weights are computed by solving a linear system
    determined by the condition that Gauss quadrature with n points
    should be exact for polynomials of degree 2n-1.

    .. cpp:function:: GaussQuadrature(unsigned int n)
    
        Create Gauss quadrature with n points

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

