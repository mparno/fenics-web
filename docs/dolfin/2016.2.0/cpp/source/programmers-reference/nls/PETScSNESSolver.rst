
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
        
    This class implements methods for solving nonlinear systems via
    PETSc's SNES interface. It includes line search and trust region
    techniques for globalising the convergence of the nonlinear
    iteration.


    .. cpp:function:: explicit PETScSNESSolver(MPI_Comm comm)
    
        Create SNES solver


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


    .. cpp:function:: void init(NonlinearProblem& nonlinear_problem, GenericVector& x)
    
        Set up the SNES object, but don't do anything yet, in case the
        user wants to access the SNES object directly


    .. cpp:function:: void set_from_options() const
    
        Set options from the PETSc options database


    .. cpp:function:: void set_options_prefix(std::string options_prefix)
    
        Sets the prefix used by PETSc when searching the PETSc options
        database


    .. cpp:function:: std::string get_options_prefix() const
    
        Returns the prefix used by PETSc when searching the PETSc
        options database


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return the MPI communicator


    .. cpp:function:: static std::vector<std::pair<std::string, std::string>> methods()
    
        Return a list of available solver methods


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: SNES snes() const
    
        Return PETSc SNES pointer


