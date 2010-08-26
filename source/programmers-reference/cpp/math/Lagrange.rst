.. Documentation for the header file dolfin/math/Lagrange.h

.. _programmers_reference_cpp_math_Mesh:

Lagrange.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Lagrange

    *Parent class*
    
        * :cpp:class:`Variable`
        
        Lagrange polynomial (basis) with given degree q determined by n = q + 1 nodal points.
        
        Example: q = 1 (n = 2)
        
        Lagrange p(1);
        p.set(0, 0.0);
        p.set(1, 1.0);
        
        This creates a Lagrange polynomial (actually two Lagrange polynomials):
        
        p(0,x) = 1 - x   (one at x = 0, zero at x = 1)
        p(1,x) = x       (zero at x = 0, one at x = 1)

    .. cpp:function:: Lagrange(const Lagrange& p)
    
        Copy constructor

    .. cpp:function:: Lagrange(unsigned int q)
    
        Constructor

    .. cpp:function:: real ddx(unsigned int i, real x)
    
        Return derivate of polynomial i at given point x

    .. cpp:function:: real dqdx(unsigned int i)
    
        Return derivative q (a constant) of polynomial

    .. cpp:function:: real eval(unsigned int i, real x)
    
        Return value of polynomial i at given point x

    .. cpp:function:: real operator() (unsigned int i, real x)
    
        Return value of polynomial i at given point x

    .. cpp:function:: real point(unsigned int i) const
    
        Return point

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: unsigned int degree() const
    
        Return degree

    .. cpp:function:: unsigned int size() const
    
        Return number of points

    .. cpp:function:: void set(unsigned int i, real x)
    
        Specify point

    .. cpp:function:: ~Lagrange()
    
        Destructor

