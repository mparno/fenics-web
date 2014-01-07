
.. Documentation for the header file dolfin/nls/PETScSNESSolver.h

.. _programmers_reference_cpp_nls_petscsnessolver:

PETScSNESSolver.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScSNESSolver

    *Parent class(es)*
    
        * :cpp:class:`PETScObject`
        
    This class implements methods for solving nonlinear systems
    via PETSc's SNES interface. It includes line search and trust
    region techniques for globalising the convergence of the
    nonlinear iteration.


    .. cpp:function:: PETScSNESSolver(std::string nls_type="default")
    
        Create SNES solver for a particular method


    .. cpp:function:: std::pair<std::size_t, bool> solve(NonlinearProblem& nonlinear_problem, GenericVector& x, const GenericVector& lb, const GenericVector& ub)
    
        Solve a nonlinear variational inequality with bound constraints
        
        *Arguments*
            nonlinear_function (:cpp:class:`NonlinearProblem`)
                The nonlinear problem.
            x (:cpp:class:`GenericVector`)
                The vector.
            lb (:cpp:class:`GenericVector`)
                The lower bound.
            ub (:cpp:class:`GenericVector`)
                The upper bound.
        
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: std::pair<std::size_t, bool> solve(NonlinearProblem& nonlinear_function, GenericVector& x)
    
        Solve abstract nonlinear problem :math:`F(x) = 0` for given
        :math:`F` and Jacobian :math:`\dfrac{\partial F}{\partial x}`.
        
        *Arguments*
            nonlinear_function (:cpp:class:`NonlinearProblem`)
                The nonlinear problem.
            x (:cpp:class:`GenericVector`)
                The vector.
        
        *Returns*
            std::pair<std::size_t, bool>
                Pair of number of Newton iterations, and whether
                iteration converged)


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: boost::shared_ptr<SNES> snes() const
    
        Return PETSc SNES pointer


    .. cpp:function:: void init(const std::string& method)
    
        Initialize SNES solver


    .. cpp:function:: void set_linear_solver_parameters()
    
        Update the linear solver parameters


    .. cpp:function:: static PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx)
    
        The callback for PETSc to compute F, the nonlinear residual


    .. cpp:function:: static PetscErrorCode FormJacobian(SNES snes, Vec x, Mat* A, Mat* B, MatStructure* flag, void* ctx)
    
        The callback for PETSc to compute A, the Jacobian


    .. cpp:function:: void set_bounds(GenericVector& x)
    
        Set the bounds on the problem from the parameters, if desired
        Here, x is passed in as a model vector from which we make our Vecs
        that tell PETSc the bounds if the "sign" parameter is used.


