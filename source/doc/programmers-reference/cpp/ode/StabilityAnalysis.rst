.. Documentation for the header file dolfin/ode/StabilityAnalysis.h

.. _programmers_reference_cpp_ode_stabilityanalysis:

StabilityAnalysis.h
===================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: StabilityAnalysis

    .. cpp:function:: StabilityAnalysis(ODE& ode, ODESolution& u)
    
        Constructor

    .. cpp:function:: void analyze_integral(uint q)
    
        Compute the integral of the q'th derivative of the dual as function of (primal) endtime T

    .. cpp:function:: void analyze_endpoint()
    
        Compute z(0) (the endpoint of the dual) as function of (primal) endtime T

