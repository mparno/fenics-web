
.. Documentation for the header file dolfin/la/PETScKrylovSolver.h

.. _programmers_reference_cpp_la_petsckrylovsolver:

PETScKrylovSolver.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
        * :cpp:class:`PETScObject`
        
    This class implements Krylov methods for linear systems of the
    form Ax = b. It is a wrapper for the Krylov solvers of PETSc.


    .. cpp:function:: enum class norm_type
    
        Norm types used in convergence testing. Not all solvers types
        support all norm types (see
        http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetNormType.html). Note
        that 'default' is a reserved keyword, so we use 'default_norm'


    .. cpp:function:: PETScKrylovSolver(MPI_Comm comm, std::string method="default", std::string preconditioner="default")
    
        Create Krylov solver for a particular method and named
        preconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method="default", std::string preconditioner="default")
    
        Create Krylov solver for a particular method and named
        preconditioner


    .. cpp:function:: PETScKrylovSolver(MPI_Comm comm, std::string method, std::shared_ptr<PETScPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: PETScKrylovSolver(std::string method, std::shared_ptr<PETScPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: PETScKrylovSolver(MPI_Comm comm, std::string method, std::shared_ptr<PETScUserPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: PETScKrylovSolver(std::string method, std::shared_ptr<PETScUserPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and
        PETScPreconditioner (shared_ptr version)


    .. cpp:function:: explicit PETScKrylovSolver(KSP ksp)
    
        Create solver wrapper of a PETSc KSP object


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(std::shared_ptr<const GenericLinearOperator> A, std::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: void set_nonzero_guess(bool nonzero_guess)
    
        Use nonzero intial guess for solution function
        (nonzero_guess=true, the solution vector x will not be zeroed
        before the solver starts)


    .. cpp:function:: void set_reuse_preconditioner(bool reuse_pc)
    
        Reuse preconditioner if true, even if matrix operator changes
        (by default preconditioner will be re-built if the matrix
        changes)


    .. cpp:function:: void set_tolerances(double relative, double absolute, double diverged, int max_iter)
    
        Set tolerances (relative residual, alsolute residial, maximum
        number of iterations)


    .. cpp:function:: void set_norm_type(norm_type type)
    
        Set norm type used in convergence testing - not all solvers
        types support all norm types


    .. cpp:function:: norm_type get_norm_type() const
    
        Get norm type used in convergence testing


    .. cpp:function:: void monitor(bool monitor_convergence)
    
        Monitor residual at each iteration


    .. cpp:function:: void set_options_prefix(std::string options_prefix)
    
        Sets the prefix used by PETSc when searching the PETSc options
        database


    .. cpp:function:: std::string get_options_prefix() const
    
        Returns the prefix used by PETSc when searching the PETSc
        options database


    .. cpp:function:: void set_from_options() const
    
        Set options from PETSc options database


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: MPI_Comm mpi_comm() const
    
        Return MPI communicator


    .. cpp:function:: KSP ksp() const
    
        Return PETSc KSP pointer


    .. cpp:function:: static std::map<std::string, std::string> methods()
    
        Return a list of available solver methods


    .. cpp:function:: static std::map<std::string, std::string> preconditioners()
    
        Return a list of available named preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::string parameter_type() const
    
        Return parameter type: "krylov_solver" or "lu_solver"


