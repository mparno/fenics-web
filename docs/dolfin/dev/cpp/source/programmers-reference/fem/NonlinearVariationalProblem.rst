
.. Documentation for the header file dolfin/fem/NonlinearVariationalProblem.h

.. _programmers_reference_cpp_fem_nonlinearvariationalproblem:

NonlinearVariationalProblem.h
=============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: NonlinearVariationalProblem

    *Parent class(es)*
    
        * :cpp:class:`Hierarchical<NonlinearVariationalProblem>`
        
    This class represents a nonlinear variational problem:
    
    Find u in V such that
    
        F(u; v) = 0  for all v in V^,
    
    where V is the trial space and V^ is the test space.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u)
    
        Create nonlinear variational problem without boundary conditions.
        The Jacobian form is not specified which requires the use of a
        nonlinear solver that does not rely on the Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const Form& J)
    
        Create nonlinear variational problem without boundary conditions.
        The Jacobian form is specified which allows the use of a nonlinear
        solver that relies on the Jacobian (using Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const BoundaryCondition& bc)
    
        Create nonlinear variational problem with a single boundary condition.
        The Jacobian form is not specified which requires the use of a
        nonlinear solver that does not rely on the Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, const BoundaryCondition& bc, const Form& J)
    
        Create nonlinear variational problem with a single boundary condition.
        The Jacobian form is specified which allows the use of a nonlinear
        solver that relies on the Jacobian (using Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, std::vector<const BoundaryCondition*> bcs)
    
        Create nonlinear variational problem with a list of boundary conditions.
        The Jacobian form is not specified which requires the use of a
        nonlinear solver that does not rely on the Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(const Form& F, Function& u, std::vector<const BoundaryCondition*> bcs, const Form& J)
    
        Create nonlinear variational problem with a list of boundary conditions.
        The Jacobian form is specified which allows the use of a nonlinear
        solver that relies on the Jacobian (using Newton's method).


    .. cpp:function:: NonlinearVariationalProblem(boost::shared_ptr<const Form> F, boost::shared_ptr<Function> u, std::vector<boost::shared_ptr<const BoundaryCondition> > bcs)
    
        Create nonlinear variational problem, shared pointer version.
        The Jacobian form is not specified which requires the use of a
        nonlinear solver that does not rely on the Jacobian.


    .. cpp:function:: NonlinearVariationalProblem(boost::shared_ptr<const Form> F, boost::shared_ptr<Function> u, std::vector<boost::shared_ptr<const BoundaryCondition> > bcs, boost::shared_ptr<const Form> J)
    
        Create nonlinear variational problem, shared pointer version.
        The Jacobian form is specified which allows the use of a nonlinear
        solver that relies on the Jacobian (using Newton's method).


    .. cpp:function:: boost::shared_ptr<const Form> residual_form() const
    
        Return residual form


    .. cpp:function:: boost::shared_ptr<const Form> jacobian_form() const
    
        Return Jacobian form


    .. cpp:function:: boost::shared_ptr<Function> solution()
    
        Return solution variable


    .. cpp:function:: boost::shared_ptr<const Function> solution() const
    
        Return solution variable (const version)


    .. cpp:function:: std::vector<boost::shared_ptr<const BoundaryCondition> > bcs() const
    
        Return boundary conditions


    .. cpp:function:: boost::shared_ptr<const FunctionSpace> trial_space() const
    
        Return trial space


    .. cpp:function:: boost::shared_ptr<const FunctionSpace> test_space() const
    
        Return test space


    .. cpp:function:: bool has_jacobian() const
    
        Check whether Jacobian has been defined


