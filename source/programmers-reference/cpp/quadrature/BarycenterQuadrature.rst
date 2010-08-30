.. Documentation for the header file dolfin/quadrature/BarycenterQuadrature.h

.. _programmers_reference_cpp_quadrature_barycenterquadrature:

BarycenterQuadrature.h
======================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: BarycenterQuadrature

    This class computes the barycenter of an arbitrary polyhedron or
    polygon in 3D and therefore allows for barycenter quadrature on
    complex polyhedrons. Note: barycenter quadrature is exact for
    polynom deg <= 1.

    .. cpp:function:: BarycenterQuadrature(const Nef_polyhedron_3& polyhedron)
    
        Create barycenter quadrature rule for given polyhedron

    .. cpp:function:: const std::vector<Point>& points() const
    
        Return points

    .. cpp:function:: const std::vector<double>& weights() const
    
        Return weights

    .. cpp:function:: uint size() const
    
        Return number of quadrature points/weights

    .. cpp:function:: void compute_quadrature(const Nef_polyhedron_3 &)
    
        Computes barycenter and weight.

