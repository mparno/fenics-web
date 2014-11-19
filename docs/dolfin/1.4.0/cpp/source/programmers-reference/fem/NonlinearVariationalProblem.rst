
.. Documentation for the header file dolfin/fem/NonlinearVariationalProblem.h

.. _programmers_reference_cpp_fem_nonlinearvariationalproblem:

NonlinearVariationalProblem.h
=============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: NonlinearVariationalProblem

    This class represents a nonlinear variational problem:
    
    Find u in V such that
    
        F(u; v) = 0  for all v in V^,
    
    where V is the trial space and V^ is the test space.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u)
    
        Create nonlinear variational problem without boundary
        conditions.  The Jacobian form is not specified which requires
        the use of a nonlinear solver that does not rely on the
        Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const Form& J)
    
        Create nonlinear variational problem without boundary
        conditions.  The Jacobian form is specified which allows the
        use of a nonlinear solver that relies on the Jacobian (using
        Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const DirichletBC& bc)
    
        Create nonlinear variational problem with a single boundary
        condition.  The Jacobian form is not specified which requires
        the use of a nonlinear solver that does not rely on the
        Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const DirichletBC& bc, const Form& J)
    
        Create nonlinear variational problem with a single boundary
        condition.  The Jacobian form is specified which allows the
        use of a nonlinear solver that relies on the Jacobian (using
        Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, std::vector<const DirichletBC*> bcs)
    
        Create nonlinear variational problem with a list of boundary
        conditions.  The Jacobian form is not specified which requires
        the use of a nonlinear solver that does not rely on the
        Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, std::vector<const DirichletBC*> bcs, const Form& J)
    
        Create nonlinear variational problem with a list of boundary
        conditions.  The Jacobian form is specified which allows the
        use of a nonlinear solver that relies on the Jacobian (using
        Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(std::shared_ptr<const Form> F, std::shared_ptr<Function> u, std::vector<std::shared_ptr<const DirichletBC> > bcs)
    
        Create nonlinear variational problem, shared pointer version.
        The Jacobian form is not specified which requires the use of a
        nonlinear solver that does not rely on the Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(std::shared_ptr<const Form> F, std::shared_ptr<Function> u, std::vector<std::shared_ptr<const DirichletBC> > bcs, std::shared_ptr<const Form> J)
    
        Create nonlinear variational problem, shared pointer version.
        The Jacobian form is specified which allows the use of a
        nonlinear solver that relies on the Jacobian (using Newton's
        method).


    .. cpp:function:: void set_bounds(std::shared_ptr<const GenericVector> lb, std::shared_ptr<const GenericVector> ub)
    
        Set the bounds for bound constrained solver


    .. cpp:function:: void set_bounds(const GenericVector& lb, const GenericVector& ub)
    
        Set the bounds for bound constrained solver


    .. cpp:function:: void set_bounds(std::shared_ptr<const Function> lb_func, std::shared_ptr<const Function> ub_func)
    
        Set the bounds for bound constrained solver


    .. cpp:function:: void set_bounds(const Function& lb_func, const Function& ub_func)
    
        Set the bounds for bound constrained solver


    .. cpp:function:: std::shared_ptr<const Form> residual_form() const
    
        Return residual form


    .. cpp:function:: std::shared_ptr<const Form> jacobian_form() const
    
        Return Jacobian form


    .. cpp:function:: std::shared_ptr<Function> solution()
    
        Return solution variable


    .. cpp:function:: std::shared_ptr<const Function> solution() const
    
        Return solution variable (const version)


    .. cpp:function:: std::vector<std::shared_ptr<const DirichletBC> > bcs() const
    
        Return boundary conditions


    .. cpp:function:: std::shared_ptr<const FunctionSpace> trial_space() const
    
        Return trial space


    .. cpp:function:: std::shared_ptr<const FunctionSpace> test_space() const
    
        Return test space


    .. cpp:function:: std::shared_ptr<const GenericVector> lower_bound() const
    
        Return lower bound


    .. cpp:function:: std::shared_ptr<const GenericVector> upper_bound() const
    
        Return upper bound


    .. cpp:function:: bool has_jacobian() const
    
        Check whether Jacobian has been defined


    .. cpp:function:: bool has_lower_bound() const
    
        Check whether lower bound has been defined


    .. cpp:function:: bool has_upper_bound() const
    
        Check whether upper bound have has defined


