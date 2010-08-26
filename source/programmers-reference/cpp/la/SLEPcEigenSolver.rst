.. Documentation for the header file dolfin/la/SLEPcEigenSolver.h

.. _programmers_reference_cpp_la_slepceigensolver:

SLEPcEigenSolver.h
==================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: class PETScMatrix
    
        Forward declarations

.. cpp:class:: SLEPcEigenSolver

    *Parent class*
    
        * :cpp:class:`Variable,`
        
        This class provides an eigenvalue solver for PETSc matrices.
        It is a wrapper for the SLEPc eigenvalue solver.
        
        The following parameters may be specified to control the solver.
        
        1. "eigenvalue spectrum"
        
        This parameter controls which part of the spectrum to compute.
        Possible values are
        
        "largest magnitude"   (eigenvalues with largest magnitude)
        "smallest magnitude"  (eigenvalues with smallest magnitude)
        "largest real"        (eigenvalues with largest double part)
        "smallest real"       (eigenvalues with smallest double part)
        "largest imaginary"   (eigenvalues with largest imaginary part)
        "smallest imaginary"  (eigenvalues with smallest imaginary part)
        "default spectrum"    (default spectrum)
        
        2. "eigenvalue solver"
        
        This parameter controls which algorithm is used by SLEPc.
        Possible values are
        
        "power"               (power iteration)
        "subspace"            (subspace iteration)
        "arnoldi"             (Arnoldi)
        "lanczos"             (Lanczos)
        "krylov-schur"        (Krylov-Schur)
        "lapack"              (LAPACK, all values, direct, only for small systems)
        "default"             (default algorithm)
        
        3. "eigenvalue tolerance"
        
        This parameter controls the tolerance used by SLEPc.
        Possible values are positive double numbers.
        
        4. "eigenvalue iterations"
        
        This parameter controls the maximum number of iterations used by SLEPc.
        Possible values are positive integers.
        
        Note that both the tolerance and the number of iterations must be
        specified if either one is specified.

    .. cpp:function:: SLEPcEigenSolver()
    
        Create eigenvalue solver

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: void get_eigenpair(double& lr, double& lc, PETScVector& r, PETScVector& c)
    
        Get the first eigenpair

    .. cpp:function:: void get_eigenpair(double& lr, double& lc, PETScVector& r, PETScVector& c, uint i)
    
        Get eigenpair i

    .. cpp:function:: void get_eigenvalue(double& lr, double& lc)
    
        Get the first eigenvalue

    .. cpp:function:: void get_eigenvalue(double& lr, double& lc, uint i)
    
        Get eigenvalue i

    .. cpp:function:: void read_parameters()
    
        Callback for changes in parameter values

    .. cpp:function:: void solve(const PETScMatrix& A)
    
        Compute all eigenpairs of the matrix A (solve Ax = \lambda x)

    .. cpp:function:: void solve(const PETScMatrix& A, const PETScMatrix& B)
    
        Compute all eigenpairs of the generalised problem Ax = \lambda Bx

    .. cpp:function:: void solve(const PETScMatrix& A, const PETScMatrix& B, uint n)
    
        Compute the n first eigenpairs of the generalised problem Ax = \lambda Bx

    .. cpp:function:: void solve(const PETScMatrix& A, uint n)
    
        Compute the n first eigenpairs of the matrix A (solve Ax = \lambda x)

    .. cpp:function:: void solve(const PETScMatrix* A, const PETScMatrix* B, uint n)
    
        Compute eigenpairs

    .. cpp:function:: ~SLEPcEigenSolver()
    
        Destructor

