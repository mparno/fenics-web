.. Documentation for the header file dolfin/ode/ODE.h

.. _programmers_reference_cpp_ode_ode:

ODE.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: ODE

    *Parent class*
    
        * :cpp:class:`Variable`
        
    An ODE represents an initial value problem of the form
    
        u'(t) = f(u(t), t) on [0, T],
    
        u(0)  = u0,
    
    where u(t) is a vector of length N.
    
    To define an ODE, a user must create a subclass of ODE and
    create the function u0() defining the initial condition, as well
    the function f() defining the right-hand side.
    
    DOLFIN provides two types of ODE solvers: a set of standard
    mono-adaptive solvers with equal adaptive time steps for all
    components as well as a set of multi-adaptive solvers with
    individual and adaptive time steps for the different
    components. The right-hand side f() is defined differently for
    the two sets of methods, with the multi-adaptive solvers
    requiring a component-wise evaluation of the right-hand
    side. Only one right-hand side function f() needs to be defined
    for use of any particular solver.
    
    It is also possible to solve implicit systems of the form
    
        M(u(t), t) u'(t) = f(u(t),t) on (0,T],
    
        u(0)  = u0,
    
    by setting the option "implicit" to true and defining the
    function M().
    
    Two different solve() functions are provided, one to solve the
    ODE on the time interval [0, T], including the solution of a
    dual problem for error control:
    
        ode.solve();
    
    Alternatively, a time interval may be given in which case the
    solution will be computed in a single sweep over the given time
    interval without solution of dual problems:
    
        ode.solve(t0, t1);
    
    This mode allows the state to be specified and retrieved in
    between intervals by calling set_state() and get_state().

    .. cpp:function:: ODE(uint N, real T)
    
        Create an ODE of size N with final time T

    .. cpp:function:: void u0(real* u) = 0
    
        Set initial values

    .. cpp:function:: void f(const real* u, real t, real* y)
    
        Evaluate right-hand side y = f(u, t), mono-adaptive version (default, optional)

    .. cpp:function:: real f(const real* u, real t, uint i)
    
        Evaluate right-hand side f_i(u, t), multi-adaptive version (optional)

    .. cpp:function:: void M(const real* dx, real* dy, const real* u, real t)
    
        Compute product dy = M dx for implicit system (optional)

    .. cpp:function:: void J(const real* dx, real* dy, const real* u, real t)
    
        Compute product dy = J dx for Jacobian J (optional)

    .. cpp:function:: void JT(const real* dx, real* dy, const real* u, real t)
    
        Compute product dy = tranpose(J) dx for Jacobian J (optional, for dual problem)

    .. cpp:function:: real dfdu(const real* u, real t, uint i, uint j)
    
        Compute entry of Jacobian (optional)

    .. cpp:function:: real timestep(real t, real k0) const
    
        Time step to use for the whole system at a given time t (optional)

    .. cpp:function:: real timestep(real t, uint i, real k0) const
    
        Time step to use for a given component at a given time t (optional)

    .. cpp:function:: bool update(const real* u, real t, bool end)
    
        Update ODE, return false to stop (optional)

    .. cpp:function:: void save(Sample& sample)
    
        Save sample (optional)

    .. cpp:function:: uint size() const
    
        Return number of components N

    .. cpp:function:: real time() const
    
        Return current time

    .. cpp:function:: real time(real t) const
    
        Return real time (might be flipped backwards for dual)

    .. cpp:function:: real endtime() const
    
        Return end time (final time T)

    .. cpp:function:: void sparse()
    
        Automatically detect sparsity (optional)

    .. cpp:function:: void solve()
    
        Solve ODE on [0, T]

    .. cpp:function:: void solve(real t0, real t1)
    
        Solve ODE on [t0, t1]

    .. cpp:function:: void solve(ODESolution& u)
    
        Solve ODE on [0, T]. Save solution in u

    .. cpp:function:: void solve(ODESolution& u, real t0, real t1)
    
        Solve ODE on [t0, t1]. Save solution in u

    .. cpp:function:: void solve_dual(ODESolution& u)
    
        Solve dual problem given an approximate solution u of the primal problem

    .. cpp:function:: void solve_dual(ODESolution& u, ODESolution& z)
    
        Solve dual and save soution in z

    .. cpp:function:: void analyze_stability(uint q, ODESolution& u)
    
        Compute stability factors as function of T (including solving the dual problem).
        The stability factor is the integral of the norm of the q'th derivative of the dual.

    .. cpp:function:: void analyze_stability_discretization(ODESolution& u)
    
        Compute stability factors as function of T (including solving the dual problem).
        The stability factor accounts for stability wrt the discretization scheme.

    .. cpp:function:: void analyze_stability_computation(ODESolution& u)
    
        Compute stability factors as function of T (including solving the dual problem).
        The stability factor accounts for stability wrt the round-off errors.

    .. cpp:function:: void analyze_stability_initial(ODESolution& u)
    
        Compute stability factors as function of T (including solving the dual problem).
        The stability factor accounts for stability wrt errors in initial data.

    .. cpp:function:: void set_state(const real* u)
    
        Set state for ODE (only available during interval stepping)

    .. cpp:function:: void get_state(real* u)
    
        Get state for ODE (only available during interval stepping)

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

