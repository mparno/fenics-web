.. Documentation for the header file dolfin/ode/ODECollection.h

.. _programmers_reference_cpp_ode_odecollection:

ODECollection.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: ODECollection

    An ODECollection represents a collection of initial value
    problems of the form
    
        u'(t) = f(u(t), t) on (0, T],
    
        u(0)  = u0,
    
    where u(t) is a vector of length N.
    
    Each ODE is governed by the same equation but a separate
    state is maintained for each ODE. Using ODECollection is
    recommended when solving a large number of ODEs and the
    overhead of instantiating a large number of ODE objects
    should be avoided.

    .. cpp:function:: ODECollection(ODE& ode, uint num_systems)
    
        Create a collection of ODEs

    .. cpp:function:: void solve(real t0, real t1)
    
        Solve ODE collection on [t0, t1]

    .. cpp:function:: void set_state(uint system, const Array<real>& u)
    
        Set state for given ODE system

    .. cpp:function:: void set_state(const Array<real>& u)
    
        Set states for all ODE systems

    .. cpp:function:: void get_state(uint system, Array<real>& u)
    
        Get state for given ODE system

    .. cpp:function:: void get_state(Array<real>& u)
    
        Get states for all ODE systems

    .. cpp:function:: void update(Array<real>& u, real t, uint system)
    
        Optional user-defined update, called between solves

