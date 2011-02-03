.. Documentation for the header file dolfin/math/basic.h

.. _programmers_reference_cpp_math_basic:

basic.h
=======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: double sqr(double x)
    
        Return the square of x

    .. cpp:function:: uint ipow(uint a, uint n)
    
        Return a to the power n

    .. cpp:function:: double rand()
    
        Return a random number, uniformly distributed between [0.0, 1.0)

    .. cpp:function:: void seed(unsigned int s)
    
        Seed random number generator

    .. cpp:function:: bool near(double x, double x0)
    
        Return true if x is within DOLFIN_EPS of x0

    .. cpp:function:: bool between(double x0, double x, double x1)
    
        Return true if x is between x0 and x1 (inclusive)

