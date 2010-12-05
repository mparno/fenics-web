.. Documentation for the header file dolfin/common/real.h

.. _programmers_reference_cpp_common_real:

real.h
======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: real real_sqrt(const real& a)
    
        Square root

    .. cpp:function:: real real_pi()
    
        Compute pi

    .. cpp:function:: void real_mat_exp(uint n, real* res, const real* A, const uint p=6)
    
        Compute matrix exponential using Pade approximation og degree p

    .. cpp:function:: real real_exp(real x)
    
        Exponential function (note: not full precision!)

    .. cpp:function:: real real_log(const real& x)
    
        Logarithmic function (note: not full precision!)

