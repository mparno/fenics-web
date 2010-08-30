.. Documentation for the header file dolfin/ode/dGqMethod.h

.. _programmers_reference_cpp_ode_dgqmethod:

dGqMethod.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: dGqMethod

    *Parent class*
    
        * :cpp:class:`Method`
        
    Contains all numeric constants, such as nodal points and nodal weights,
    needed for the dG(q) method. The order q must be at least 0. Note that
    q refers to the polynomial order and not the order of convergence for
    the method, which is 2q + 1.

    .. cpp:function:: real error(real k, real r) const
    
        Compute error estimate (modulo stability factor)

    .. cpp:function:: real residual(real x0, real values[], real f, real k) const
    
        Compute residual at right end-point

    .. cpp:function:: real timestep(real r, real tol, real k0, real kmax) const
    
        Compute new time step based on the given residual

    .. cpp:function:: real ueval(real x0, real values[], real tau) const
    
        Evaluate solution at given point

    .. cpp:function:: real ueval(real x0, real values[], uint i) const
    
        Evaluate solution at given node (inline optimized)

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

