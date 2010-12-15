.. Documentation for the header file dolfin/ode/ODESolution.h

.. _programmers_reference_cpp_ode_odesolution:

ODESolution.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: ODESolutionData

.. cpp:class:: ODESolution

    .. cpp:function:: ODESolution()
    
        Create solution data for given ODE

    .. cpp:function:: void flush()
    
        Make object ready for evaluating, set to read mode

    .. cpp:function:: void eval(const real& t, real* y)
    
        Evaluate (interpolate) value of solution at given time

    .. cpp:function:: ODESolutionData& get_timeslab(uint index)
    
        Get timeslab (used when iterating)

    .. cpp:function:: const real* get_weights() const
    
        Get pointer to weights

.. cpp:class:: ODESolutionIterator

