.. Documentation for the header file dolfin/ode/StabilityAnalysis.h

.. _programmers_reference_cpp_ode_stabilityanalysis:

StabilityAnalysis.h
===================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: class ODESolution
    
        This class computes the stabilityfactors as a function of time
        S(t). As the stabilityfactors are defined it should solve the dual
        for each timestep. However, we can take advantage of the dual
        being linear.

.. cpp:class:: StabilityAnalysis

    .. cpp:function:: StabilityAnalysis(ODE& ode, ODESolution& u)
    
        Constructor

    .. cpp:function:: void analyze_endpoint()
    
        Compute z(0) (the endpoint of the dual) as function of (primal) endtime T

    .. cpp:function:: void analyze_integral(uint q)
    
        Compute the integral of the q'th derivative of the dual as function of (primal) endtime T

    .. cpp:function:: ~StabilityAnalysis()
    
        Destructor

