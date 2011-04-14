.. Documentation for the header file dolfin/quadrature/GaussianQuadrature.h

.. _programmers_reference_cpp_quadrature_gaussianquadrature:

GaussianQuadrature.h
====================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GaussianQuadrature

    *Parent class*
    
        * :cpp:class:`Quadrature`
        
    Gaussian-type quadrature rule on the double line,
    including Gauss, Radau, and Lobatto quadrature.
    
    Points and weights are computed to be exact within a tolerance
    of DOLFIN_EPS. Comparing with known exact values for n <= 3 shows
    that we obtain full precision (16 digits, error less than 2e-16).

