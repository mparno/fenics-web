
.. Documentation for the header file dolfin/math/Lagrange.h

.. _programmers_reference_cpp_math_lagrange:

Lagrange.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Lagrange

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    Lagrange polynomial (basis) with given degree q determined by
    n = q + 1 nodal points.
    
    Example: q = 1 (n = 2)
    
      Lagrange p(1);
      p.set(0, 0.0);
      p.set(1, 1.0);
    
    It is the callers reponsibility that the points are distinct.
    
    This creates a Lagrange polynomial (actually two Lagrange
    polynomials):
    
      p(0,x) = 1 - x   (one at x = 0, zero at x = 1)
      p(1,x) = x       (zero at x = 0, one at x = 1)
    


    .. cpp:function:: Lagrange(std::size_t q)
    
        Constructor


    .. cpp:function:: Lagrange(const Lagrange& p)
    
        Copy constructor


    .. cpp:function:: void set(std::size_t i, double x)
    
        Specify point


    .. cpp:function:: std::size_t size() const
    
        Return number of points


    .. cpp:function:: std::size_t degree() const
    
        Return degree


    .. cpp:function:: double point(std::size_t i) const
    
        Return point


    .. cpp:function:: double operator() (std::size_t i, double x)
    
        Return value of polynomial i at given point x


    .. cpp:function:: double eval(std::size_t i, double x)
    
        Return value of polynomial i at given point x


    .. cpp:function:: double ddx(std::size_t i, double x)
    
        Return derivate of polynomial i at given point x


    .. cpp:function:: double dqdx(std::size_t i)
    
        Return derivative q (a constant) of polynomial


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


