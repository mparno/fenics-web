.. Documentation for the header file dolfin/ode/ComplexODE.h

.. _programmers_reference_cpp_ode_Mesh:

ComplexODE.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: ComplexODE

    *Parent class*
    
        * :cpp:class:`ODE`
        
        A ComplexODE represents an initial value problem for a system of
        complex-valued ordinary differential equations:
        
        M(z, t) z'(t) = f(z(t), t) on (0,T]
        
        z(0) = z0,
        
        where z(t) is a complex-valued vector of length n. The imaginary
        unit is provided by the member variable j satisfying j^2 = -1.
        
        This class is a wrapper for a standard real-valued ODE, and
        provides an interface that automatically translates the given
        complex-valued ODE of size n to a standard real-valued ODE of
        size N = 2n.
        
        The double and imaginary parts of the solution are stored in the
        following order in the solution vector u(t):
        
        u = (Re z0, Im z0, Re z1, Im z1, ..., Re z_n-1, Im z_n-1).

    .. cpp:function:: ComplexODE(uint n, real T)
    
        Constructor

    .. cpp:function:: bool update(const real* u, real t, bool end)
    
        Update for real-valued ODE

    .. cpp:function:: real f(const real* u, real t, uint i)
    
        Return right-hand side for real-valued ODE

    .. cpp:function:: real timestep(uint i)
    
        Return time step for real-valued ODE

    .. cpp:function:: virtual bool update(const complex z[], real t, bool end)
    
        Update ODE, return false to stop (optional)

    .. cpp:function:: virtual complex f(const complex z[], real t, uint i)
    
        Evaluate right-hand side (multi-adaptive version)

    .. cpp:function:: virtual real k(uint i)
    
        Return time step for component i (optional)

    .. cpp:function:: virtual void J(const complex x[], complex y[], const complex u[], real t)
    
        Compute product y = Jx for Jacobian J

    .. cpp:function:: virtual void M(const complex x[], complex y[], const complex z[], real t)
    
        Compute product y = Mx for implicit system

    .. cpp:function:: virtual void f(const complex z[], real t, complex y[])
    
        Evaluate right-hand side (mono-adaptive version)

    .. cpp:function:: virtual void z0(complex z[]) = 0
    
        Set initial values

    .. cpp:function:: void J(const real* x, real* y, const real* u, real t)
    
        Compute product y = Jx for real-valued ODE

    .. cpp:function:: void M(const real* x, real* y, const real* u, real t)
    
        Compute product y = Mx for real-valued ODE

    .. cpp:function:: void f(const real* u, real t, real* y)
    
        Evaluate right-hand side for real-valued ODE

    .. cpp:function:: void u0(real* u)
    
        Return initial value for real-valued ODE

    .. cpp:function:: ~ComplexODE()
    
        Destructor

.. cpp:class:: DummyComplexODE

    *Parent class*
    
        * :cpp:class:`ODE`
        
        Dummy implementation of ComplexODE used when DOLFIN is compiled
        with GMP support in which case ComplexODE is not available

