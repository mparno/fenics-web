
.. Documentation for the header file dolfin/la/uBLASKrylovSolver.h

.. _programmers_reference_cpp_la_ublaskrylovsolver:

uBLASKrylovSolver.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class implements Krylov methods for linear systems
    of the form Ax = b using uBLAS data types.


    .. cpp:function:: uBLASKrylovSolver(std::string method="default", std::string preconditioner="default")
    
        Create Krylov solver for a particular method and preconditioner


    .. cpp:function:: uBLASKrylovSolver(uBLASPreconditioner& pc)
    
        Create Krylov solver for a particular uBLASPreconditioner


    .. cpp:function:: uBLASKrylovSolver(std::string method, uBLASPreconditioner& pc)
    
        Create Krylov solver for a particular method and uBLASPreconditioner


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Solve the operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericLinearOperator> A, const boost::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: const GenericLinearOperator& get_operator() const
    
        Return the operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::size_t solve_krylov(const MatA& A, uBLASVector& x, const uBLASVector& b, const MatP& P)
    
        Select solver and solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solveCG(const Mat& A, uBLASVector& x, const uBLASVector& b, bool& converged) const
    
        Solve linear system Ax = b using CG


    .. cpp:function:: std::size_t solveGMRES(const Mat& A, uBLASVector& x, const uBLASVector& b, bool& converged) const
    
        Solve linear system Ax = b using restarted GMRES


    .. cpp:function:: std::size_t solveBiCGStab(const Mat& A, uBLASVector& x, const uBLASVector& b, bool& converged) const
    
        Solve linear system Ax = b using BiCGStab


    .. cpp:function:: void select_preconditioner(std::string preconditioner)
    
        Select and create named preconditioner


    .. cpp:function:: void read_parameters()
    
        Read solver parameters


