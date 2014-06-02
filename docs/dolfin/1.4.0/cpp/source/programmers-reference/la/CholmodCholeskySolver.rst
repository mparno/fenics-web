
.. Documentation for the header file dolfin/la/CholmodCholeskySolver.h

.. _programmers_reference_cpp_la_cholmodcholeskysolver:

CholmodCholeskySolver.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CholmodCholeskySolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class implements the direct solution (Cholesky
    factorization) of linear systems of the form Ax = b. Sparse
    matrices are solved using CHOLMOD
    http://www.cise.ufl.edu/research/sparse/cholmod/ if installed.


    .. cpp:function:: CholmodCholeskySolver()
    
        Constructor


    .. cpp:function:: CholmodCholeskySolver(std::shared_ptr<const GenericLinearOperator> A)
    
        Constructor


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Solve the operator (matrix)


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b for a sparse matrix using CHOLMOD


    .. cpp:function:: std::size_t factorize(const GenericLinearOperator& A)
    
        Cholesky-factor sparse matrix A if CHOLMOD is installed


    .. cpp:function:: std::size_t factorized_solve(GenericVector& x, const GenericVector& b)
    
        Solve factorized system (CHOLMOD).


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


.. cpp:class:: Cholmod

    .. cpp:function:: void clear()
    
        Clear data


    .. cpp:function:: void init(long int* Ap, long int* Ai, double* Ax, std::size_t M, std::size_t nz)
    
        Initialise with matrix


    .. cpp:function:: void factorize()
    
        Factorize


    .. cpp:function:: void factorized_solve(double*x, const double* b)
    
        Factorized solve


    .. cpp:function:: cholmod_dense* residual(cholmod_dense* x, cholmod_dense* b)
    
        Compute residual: b-Ax


    .. cpp:function:: double residual_norm(cholmod_dense* r, cholmod_dense* x, cholmod_dense* b)
    
        Compute residual norm


    .. cpp:function:: void refine_once(cholmod_dense* x, cholmod_dense* r)
    
        Perform one refinement


    .. cpp:function:: void check_status(std::string function)
    
        Check status flag returned by an CHOLMOD function


