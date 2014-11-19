
.. Documentation for the header file dolfin/fem/PointSource.h

.. _programmers_reference_cpp_fem_pointsource:

PointSource.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PointSource

    This class provides an easy mechanism for adding a point source
    (Dirac delta function) to the right-hand side vector in a
    variational problem. The associated function space must be
    scalar in order for the inner product with the (scalar) Dirac
    delta function to be well defined.


    .. cpp:function:: PointSource(const FunctionSpace& V, const Point& p, double magnitude=1.0)
    
        Create point source at given point of given magnitude


    .. cpp:function:: PointSource(boost::shared_ptr<const FunctionSpace> V, const Point& p, double magnitude=1.0)
    
        Create point source at given point of given magnitude


    .. cpp:function:: void apply(GenericVector& b)
    
        Apply (add) point source to right-hand side vector


