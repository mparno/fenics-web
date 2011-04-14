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
        
    A :cpp:class:`VariationalProblem` represents a (system of) partial
    differential equation(s) in variational form:
    
    Find u_h in V_h such that
    
        F(u_h; v) = 0  for all v in V_h',
    
    where V_h is the trial space and V_h' is the test space.
    
    The variational problem is specified in terms of a pair of
    _Form_s and, optionally, a set of _BoundaryCondition_s.
    
    The pair of forms may either specify a nonlinear problem
    
       (1) F(u_h; v) = 0
    
    in terms of the residual F and its derivative J = F':
    
       F, J  (F linear, J bilinear)
    
    or a linear problem
    
       (2) F(u_h; v) = a(u_h, v) - L(v) = 0
    
    in terms of the bilinear form a and a linear form L:
    
       a, L  (a bilinear, L linear)
    
    Thus, a pair of forms is interpreted either as a nonlinear
    problem or a linear problem depending on the ranks of the given
    forms.

    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1)
    
        Define variational problem with natural boundary conditions

    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1, const BoundaryCondition& bc)
    
        Define variational problem with a single Dirichlet boundary condition

    .. cpp:function:: VariationalProblem(const Form& form_0, const Form& form_1, const std::vector<const BoundaryCondition*>& bcs)
    
        Define variational problem with a list of Dirichlet boundary conditions

    .. cpp:function:: VariationalProblem(boost::shared_ptr<const Form> form_0, boost::shared_ptr<const Form> form_1, std::vector<boost::shared_ptr<const BoundaryCondition> > bcs)
    
        Define variational problem with a list of Dirichlet boundary conditions

    .. cpp:function:: void solve(Function& u) const
    
        Solve variational problem

    .. cpp:function:: void solve(Function& u0, Function& u1) const
    
        Solve variational problem and extract sub functions

    .. cpp:function:: void solve(Function& u0, Function& u1, Function& u2) const
    
        Solve variational problem and extract sub functions

    .. cpp:function:: void solve(Function& u, const double tol, GoalFunctional& M) const
    
        Solve variational problem adaptively to within given tolerance

    .. cpp:function:: void solve(Function& u, const double tol, Form& M, ErrorControl& ec) const
    
        Solve variational problem adaptively to within given tolerance

    .. cpp:function:: const bool is_nonlinear() const
    
        Return true if problem is non-linear

    .. cpp:function:: const FunctionSpace& trial_space() const
    
        Return trial space for variational problem

    .. cpp:function:: const FunctionSpace& test_space() const
    
        Return test space for variational problem

    .. cpp:function:: const Form& bilinear_form() const
    
        Return the bilinear form

    .. cpp:function:: boost::shared_ptr<const Form> bilinear_form_shared_ptr() const
    
        Return the bilinear form (shared_ptr version)

    .. cpp:function:: boost::shared_ptr<const Form> form_0_shared_ptr() const
    
        Return form_0 (shared_ptr version)

    .. cpp:function:: boost::shared_ptr<const Form> form_1_shared_ptr() const
    
        Return form_1 (shared_ptr version)

    .. cpp:function:: const Form& linear_form() const
    
        Return the linear form

    .. cpp:function:: boost::shared_ptr<const Form> linear_form_shared_ptr() const
    
        Return the linear form (shared_ptr version)

    .. cpp:function:: const std::vector<const BoundaryCondition*> bcs() const
    
        Return the list of boundary conditions

    .. cpp:function:: const std::vector<boost::shared_ptr<const BoundaryCondition> > bcs_shared_ptr() const
    
        Return the list of boundary conditions (shared_ptr version)

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

