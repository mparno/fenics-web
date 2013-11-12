
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


    .. cpp:function:: EpetraKrylovSolver(std::string method = "default", std::string preconditioner = "default")
    
        Create Krylov solver for a particular method and preconditioner


    .. cpp:function:: EpetraKrylovSolver(std::string method, TrilinosPreconditioner& preconditioner)
    
        Create Krylov solver for a particular method and
        TrilinosPreconditioner


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Set the operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericLinearOperator> A, const boost::shared_ptr<const GenericLinearOperator> P)
    
        Set the operator (matrix)


    .. cpp:function:: const GenericLinearOperator& get_operator() const
    
        Get the operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const EpetraMatrix& A, EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: double residual(const std::string residualtype) const
    
        Return residual from most recent solve


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<AztecOO> aztecoo() const
    
        Return pointer to Aztec00


