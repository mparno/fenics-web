
.. Documentation for the header file dolfin/la/EpetraKrylovSolver.h

.. _programmers_reference_cpp_la_epetrakrylovsolver:

EpetraKrylovSolver.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EpetraKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class implements Krylov methods for linear systems
    of the form Ax = b. It is a wrapper for the Krylov solvers
    of Epetra.


    .. cpp:function:: EpetraKrylovSolver(std::string method = "default", std::string pc_type = "default")
    
        Create Krylov solver for a particular method and preconditioner


    .. cpp:function:: EpetraKrylovSolver(std::string method, TrilinosPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and TrilinosPreconditioner


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericMatrix> A)
    
        Set the operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericMatrix> A, const boost::shared_ptr<const GenericMatrix> P)
    
        Set the operator (matrix)


    .. cpp:function:: const GenericMatrix& get_operator() const
    
        Get the operator (matrix)


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(const EpetraMatrix& A, EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<AztecOO> aztecoo() const
    
        Return pointer to Aztec00


