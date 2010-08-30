.. Documentation for the header file dolfin/quadrature/RadauQuadrature.h

.. _programmers_reference_cpp_quadrature_radauquadrature:

RadauQuadrature.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: RadauQuadrature

    *Parent class*
    
        * :cpp:class:`GaussianQuadrature`
        
    Radau (Gauss-Radau) quadrature on the interval [-1,1].
    The n quadrature points are given by the zeros of
    
        ( Pn-1(x) + Pn(x) ) / (1+x)
    
    where Pn is the n:th Legendre polynomial.
    
    The quadrature points are computed using Newton's method, and
    the quadrature weights are computed by solving a linear system
    determined by the condition that Radau quadrature with n points
    should be exact for polynomials of degree 2n-2.

    .. cpp:function:: RadauQuadrature(unsigned int n)
    
        Create Radau quadrature with n points

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

