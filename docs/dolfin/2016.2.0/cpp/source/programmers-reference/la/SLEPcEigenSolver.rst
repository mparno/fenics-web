
.. Documentation for the header file dolfin/la/SLEPcEigenSolver.h

.. _programmers_reference_cpp_la_slepceigensolver:

SLEPcEigenSolver.h
==================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SLEPcEigenSolver

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
        * :cpp:class:`PETScObject`
        
    This class provides an eigenvalue solver for PETSc matrices. It
    is a wrapper for the SLEPc eigenvalue solver.
    
    The following parameters may be specified to control the solver.
    
    1. "spectrum"
    
    This parameter controls which part of the spectrum to compute.
    Possible values are
    
      "largest magnitude"   (eigenvalues with largest magnitude)
      "smallest magnitude"  (eigenvalues with smallest magnitude)
      "largest real"        (eigenvalues with largest double part)
      "smallest real"       (eigenvalues with smallest double part)
      "largest imaginary"   (eigenvalues with largest imaginary part)
      "smallest imaginary"  (eigenvalues with smallest imaginary part)
      "target magnitude"    (eigenvalues closest to target in magnitude)
      "target real"         (eigenvalues closest to target in real part)
      "target imaginary"    (eigenvalues closest to target in imaginary part)
    
    2. "solver"
    
    This parameter controls which algorithm is used by SLEPc.
    Possible values are
    
      "power"               (power iteration)
      "subspace"            (subspace iteration)
      "arnoldi"             (Arnoldi)
      "lanczos"             (Lanczos)
      "krylov-schur"        (Krylov-Schur)
      "lapack"              (LAPACK, all values, direct, small systems only)
      "arpack"              (ARPACK)
    
    3. "tolerance"
    
    This parameter controls the tolerance used by SLEPc.  Possible
    values are positive double numbers.
    
    4. "maximum_iterations"
    
    This parameter controls the maximum number of iterations used by SLEPc.
    Possible values are positive integers.
    
    Note that both the tolerance and the number of iterations must be
    specified if either one is specified.
    
    5. "problem_type"
    
    This parameter can be used to give extra information about the
    type of the eigenvalue problem. Some solver types require this
    extra piece of information. Possible values are:
    
      "hermitian"               (Hermitian)
      "non_hermitian"           (Non-Hermitian)
      "gen_hermitian"           (Generalized Hermitian)
      "gen_non_hermitian"       (Generalized Non-Hermitian)
      "pos_gen_non_hermitian"   (Generalized Non-Hermitian with positive semidefinite B)
    
    6. "spectral_transform"
    
    This parameter controls the application of a spectral
    transform. A spectral transform can be used to enhance the
    convergence of the eigensolver and in particular to only compute
    eigenvalues in the interior of the spectrum. Possible values
    are:
    
      "shift-and-invert"      (A shift-and-invert transform)
    
    Note that if a spectral transform is given, then also a non-zero
    spectral shift parameter has to be provided.
    
    The default is no spectral transform.
    
    7. "spectral_shift"
    
    This parameter controls the spectral shift used by the spectral
    transform and must be provided if a spectral transform is
    given. The possible values are real numbers.


    .. cpp:function:: explicit SLEPcEigenSolver(MPI_Comm comm)
    
        Create eigenvalue solver


    .. cpp:function:: explicit SLEPcEigenSolver(EPS eps)
    
        Create eigenvalue solver from EPS object


    .. cpp:function:: explicit SLEPcEigenSolver(std::shared_ptr<const PETScMatrix> A)
    
        Create eigenvalue solver for Ax = \lambda


    .. cpp:function:: SLEPcEigenSolver(MPI_Comm comm, std::shared_ptr<const PETScMatrix> A)
    
        Create eigenvalue solver for Ax = \lambda x


    .. cpp:function:: SLEPcEigenSolver(std::shared_ptr<const PETScMatrix> A, std::shared_ptr<const PETScMatrix> B)
    
        Create eigenvalue solver for Ax = \lambda x on MPI_COMM_WORLD


    .. cpp:function:: SLEPcEigenSolver(MPI_Comm comm, std::shared_ptr<const PETScMatrix> A, std::shared_ptr<const PETScMatrix> B)
    
        Create eigenvalue solver for Ax = \lambda x


    .. cpp:function:: void set_operators(std::shared_ptr<const PETScMatrix> A, std::shared_ptr<const PETScMatrix> B)
    
        Set opeartors (B may be nullptr for regular eigenvalues
        problems)


    .. cpp:function:: void solve()
    
        Compute all eigenpairs of the matrix A (solve Ax = \lambda x)


    .. cpp:function:: void solve(std::size_t n)
    
        Compute the n first eigenpairs of the matrix A (solve Ax = \lambda x)


    .. cpp:function:: void get_eigenvalue(double& lr, double& lc) const
    
        Get the first eigenvalue


    .. cpp:function:: void get_eigenpair(double& lr, double& lc, GenericVector& r, GenericVector& c) const
    
        Get the first eigenpair


    .. cpp:function:: void get_eigenpair(double& lr, double& lc, PETScVector& r, PETScVector& c) const
    
        Get the first eigenpair


    .. cpp:function:: void get_eigenvalue(double& lr, double& lc, std::size_t i) const
    
        Get eigenvalue i


    .. cpp:function:: void get_eigenpair(double& lr, double& lc, GenericVector& r, GenericVector& c, std::size_t i) const
    
        Get eigenpair i


    .. cpp:function:: void get_eigenpair(double& lr, double& lc, PETScVector& r, PETScVector& c, std::size_t i) const
    
        Get eigenpair i


    .. cpp:function:: std::size_t get_iteration_number() const
    
        Get the number of iterations used by the solver


    .. cpp:function:: std::size_t get_number_converged() const
    
        Get the number of converged eigenvalues


    .. cpp:function:: void set_deflation_space(const PETScVector& deflation_space)
    
        Set deflation space


    .. cpp:function:: void set_options_prefix(std::string options_prefix)
    
        Sets the prefix used by PETSc when searching the PETSc options
        database


    .. cpp:function:: std::string get_options_prefix() const
    
        Returns the prefix used by PETSc when searching the PETSc
        options database


    .. cpp:function:: EPS eps() const
    
        Return SLEPc EPS pointer


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void read_parameters()
    
        Callback for changes in parameter values


