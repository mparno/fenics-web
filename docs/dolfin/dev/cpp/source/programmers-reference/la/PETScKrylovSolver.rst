
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
        
    This class implements Krylov methods for linear systems
    of the form Ax = b. It is a wrapper for the Krylov solvers
    of PETSc.


    .. cpp:function:: PETScKrylovSolver(std::string method = "default", std::string pc_type = "default")
    
        Create Krylov solver for a particular method and names preconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method, PETScPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and PETScPreconditioner


    .. cpp:function:: PETScKrylovSolver(std::string method, PETScUserPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and PETScPreconditioner


    .. cpp:function:: explicit PETScKrylovSolver(boost::shared_ptr<KSP> ksp)
    
        Create solver from given PETSc KSP pointer


    .. cpp:function:: void set_operator(const GenericMatrix& A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operator(const PETScBaseMatrix& A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const GenericMatrix& A, const GenericMatrix& P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: void set_operators(const PETScBaseMatrix& A, const PETScBaseMatrix& P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: const GenericMatrix& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(const PETScBaseMatrix& A, PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<KSP> ksp() const
    
        Return PETSc KSP pointer


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: void init(const std::string& method)
    
        Initialize KSP solver


    .. cpp:function:: void write_report(int num_iterations, KSPConvergedReason reason)
    
        Report the number of iterations


