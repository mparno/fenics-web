.. Documentation for the header file dolfin/ode/ODESolution.h

.. _programmers_reference_cpp_ode_odesolution:

ODESolution.h
=============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: class Lagrange
    
        ODESolution stores the solution from the ODE solver, primarily to
        be able to solve the dual problem. A number of interpolated values
        is cached, since the ODE solver repeatedly requests evaluation of
        the same t.
        
        The samples are stored in memory if possible, otherwise stored
        in a temporary file and fetched from disk in blocks when needed.
        
        Since GMP at the moment doesn't support saving binary operands
        on disk this class uses ascii files for storage.
        Fortunately storing operands on disk in binary is planned in
        the next major release of GMP.

.. cpp:class:: ODESolutionData

.. cpp:class:: ODESolution

    .. cpp:function:: ODESolution()
    
        Create solution data for given ODE

    .. cpp:function:: ODESolutionData& get_timeslab(uint index)
    
        Get timeslab (used when iterating)

    .. cpp:function:: const real* get_weights() const
    
        Get pointer to weights

    .. cpp:function:: void eval(const real& t, real* y)
    
        Evaluate (interpolate) value of solution at given time

    .. cpp:function:: void flush()
    
        Make object ready for evaluating, set to read mode

.. cpp:class:: ODESolutionIterator

