.. Documentation for the header file dolfin/quadrature/LobattoQuadrature.h

.. _programmers_reference_cpp_quadrature_Mesh:

LobattoQuadrature.h
===================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: LobattoQuadrature

    *Parent class*
    
        * :cpp:class:`GaussianQuadrature`
        
        Lobatto (Gauss-Lobatto) quadrature on the interval [-1,1].
        The n quadrature points are given by the end-points -1 and 1,
        and the zeros of P{n-1}'(x), where P{n-1}(x) is the (n-1):th
        Legendre polynomial.
        
        The quadrature points are computed using Newton's method, and
        the quadrature weights are computed by solving a linear system
        determined by the condition that Lobatto quadrature with n points
        should be exact for polynomials of degree 2n-3.

    .. cpp:function:: LobattoQuadrature(unsigned int n)
    
        Create Lobatto quadrature with n points

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

