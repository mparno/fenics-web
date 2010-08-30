.. Documentation for the header file dolfin/fem/VariationalProblem.h

.. _programmers_reference_cpp_fem_variationalproblem:

VariationalProblem.h
====================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: VariationalProblem

    *Parent class*
    
        * :cpp:class:`Variable,`
        
    This class represents a (system of) partial differential
    equation(s) in variational form: Find u in V such that
    
        F_u(v) = 0  for all v in V'.
    
    The variational problem is defined in terms of a bilinear
    form a(v, u) and a linear for L(v).
    
    For a linear variational problem, F_u(v) = a(v, u) - L(v),
    the forms should correspond to the canonical formulation
    
        a(v, u) = L(v)  for all v in V'.
    
    For a nonlinear variational problem, the forms should
    be given by
    
        a(v, u) = F_u'(v) u = F_u'(v, u),
        L(v)    = F(v),
    
    that is, a(v, u) should be the Frechet derivative of F_u
    with respect to u, and L = F.
    
    Parameters:
    
        "linear solvers": "direct" or "iterative" (default: "direct")
        "symmetric":      true or false (default: false)

    .. cpp:function:: NewtonSolver& newton_solver()
    
        Return Newton solver (only useful when solving a nonlinear problem)

    .. cpp:function:: VariationalProblem(const Form& a,
                                         const Form& L,
                                         const BoundaryCondition& bc,
                                         bool nonlinear=false)
    
        Define variational problem with a single Dirichlet boundary conditions

    .. cpp:function:: VariationalProblem(const Form& a,
                                         const Form& L,
                                         const std::vector<const BoundaryCondition*>& bcs,
                                         bool nonlinear=false)
    
        Define variational problem with a list of Dirichlet boundary conditions

    .. cpp:function:: VariationalProblem(const Form& a,
                                         const Form& L,
                                         const std::vector<const BoundaryCondition*>& bcs,
                                         const MeshFunction<uint>* cell_domains,
                                         const MeshFunction<uint>* exterior_facet_domains,
                                         const MeshFunction<uint>* interior_facet_domains,
                                         bool nonlinear=false)
    
        Define variational problem with a list of Dirichlet boundary conditions
        and subdomains

    .. cpp:function:: VariationalProblem(const Form& a,
                       const Form& L,
                       bool nonlinear=false)
    
        Define variational problem with natural boundary conditions

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: void F(GenericVector& b, const GenericVector& x)
    
        Compute F at current point x

    .. cpp:function:: void J(GenericMatrix& A, const GenericVector& x)
    
        Compute J = F' at current point x

    .. cpp:function:: void solve(Function& u)
    
        Solve variational problem

    .. cpp:function:: void solve(Function& u0, Function& u1)
    
        Solve variational problem and extract sub functions

    .. cpp:function:: void solve(Function& u0, Function& u1, Function& u2)
    
        Solve variational problem and extract sub functions

    .. cpp:function:: void update(const GenericVector& x)
    
        Optional callback called before calls to F() and J()

