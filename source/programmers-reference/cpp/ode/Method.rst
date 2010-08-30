.. Documentation for the header file dolfin/ode/Method.h

.. _programmers_reference_cpp_ode_method:

Method.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Method

    *Parent class*
    
        * :cpp:class:`Variable`
        
    Base class for cGqMethod and dGqMethod, which contain all numeric constants,
    such as nodal points and nodal weights, needed for the method.

    .. cpp:function:: Method(unsigned int q, unsigned int nq, unsigned int nn)
    
        Constructor

    .. cpp:function:: Type type() const
    
        Return type (inline optimized)

    .. cpp:function:: const Lagrange get_trial() const
    
        Get trial functions

    .. cpp:function:: const real* get_quadrature_weights() const
    
        Get quadrature weights

    .. cpp:function:: real derivative(unsigned int i) const
    
        Evaluation of derivative of basis function i at t = 1 (inline optimized)

    .. cpp:function:: real error(real k, real r) const = 0
    
        Compute error estimate (modulo stability factor)

    .. cpp:function:: real eval(unsigned int i, real tau) const
    
        Evaluation of trial space basis function i at given tau (inline optimized)

    .. cpp:function:: real npoint(unsigned int i) const
    
        Return nodal point (inline optimized)

    .. cpp:function:: real nweight(unsigned int i, unsigned int j) const
    
        Return nodal weight j for node i, including quadrature (inline optimized)

    .. cpp:function:: real qpoint(unsigned int i) const
    
        Return quadrature point (inline optimized)

    .. cpp:function:: real qweight(unsigned int i) const
    
        Return quadrature weight, including only quadrature (inline optimized)

    .. cpp:function:: real residual(real x0, real values[], real f, real k) const = 0
    
        Compute residual at right end-point

    .. cpp:function:: real timestep(real r, real tol, real k0, real kmax) const = 0
    
        Compute new time step based on the given residual

    .. cpp:function:: real ueval(real x0, real values[], real tau) const = 0
    
        Evaluate solution at given point

    .. cpp:function:: real ueval(real x0, real values[], uint i) const = 0
    
        Evaluate solution at given node

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: unsigned int degree() const
    
        Return degree (inline optimized)

    .. cpp:function:: unsigned int nsize() const
    
        Return number of nodal points (inline optimized)

    .. cpp:function:: unsigned int order() const
    
        Return order (inline optimized)

    .. cpp:function:: unsigned int qsize() const
    
        Return number of quadrature points (inline optimized)

    .. cpp:function:: void get_nodal_values(const real& x0, const real* x, real* nodal_values) const = 0
    
        Get nodal values

    .. cpp:function:: void update(real x0, real f[], real k, real values[]) const
    
        Update solution values using fixed-point iteration

    .. cpp:function:: void update(real x0, real f[], real k, real values[], real alpha) const
    
        Update solution values using fixed-point iteration (damped version)

