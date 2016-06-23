
.. Documentation for the header file dolfin/math/basic.h

.. _programmers_reference_cpp_math_basic:

basic.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: std::size_t ipow(std::size_t a, std::size_t n)

    Return a to the power n.
    NOTE: Overflow is not checked!


.. cpp:function:: double rand()

    Return a random number, uniformly distributed between [0.0, 1.0)


.. cpp:function:: void seed(std::size_t s)

    Seed random number generator


.. cpp:function:: bool near(double x, double x0, double eps=DOLFIN_EPS)

    Check whether x is close to x0 (to within DOLFIN_EPS)


.. cpp:function:: bool between(double x, std::pair<double, double> range)

    Check whether x is between x0 and x1 (inclusive, to within DOLFIN_EPS)


