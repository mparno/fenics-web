.. Documentation for the header file dolfin/ode/cGqMethod.h

.. _programmers_reference_cpp_ode_Mesh:

cGqMethod.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: cGqMethod

    *Parent class*
    
        * :cpp:class:`Method`
        
        Contains all numeric constants, such as nodal points and nodal weights,
        needed for the cG(q) method. The order q must be at least 1. Note that
        q refers to the polynomial order and not the order of convergence for
        the method, which is 2q.

    .. cpp:function:: inline real ueval(real x0, real values[], uint i) const
    
        Evaluate solution at given node (inline optimized)

    .. cpp:function:: real error(real k, real r) const
    
        Compute error estimate (modulo stability factor)

    .. cpp:function:: real residual(real x0, real values[], real f, real k) const
    
        Compute residual at right end-point

    .. cpp:function:: real timestep(real r, real tol, real k0, real kmax) const
    
        Compute new time step based on the given residual

    .. cpp:function:: real ueval(real x0, real values[], real tau) const
    
        Evaluate solution at given point

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: void get_nodal_values(const real& x0, const real* x, real* nodal_values) const
    
        Replace the solution values with the nodal values solution polynomial

