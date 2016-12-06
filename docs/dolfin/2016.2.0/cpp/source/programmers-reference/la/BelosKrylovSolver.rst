
.. Documentation for the header file dolfin/la/BelosKrylovSolver.h

.. _programmers_reference_cpp_la_beloskrylovsolver:

BelosKrylovSolver.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BelosKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class implements Krylov methods for linear systems
    of the form Ax = b. It is a wrapper for the Belos iterative solver
    from Trilinos


    .. cpp:function:: BelosKrylovSolver(std::string method = "default", std::string preconditioner = "default")
    
        Create Krylov solver for a particular method and names
        preconditioner


    .. cpp:function:: BelosKrylovSolver(std::string method, std::shared_ptr<TrilinosPreconditioner> preconditioner)
    
        Create Krylov solver for a particular method and TrilinosPreconditioner


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(std::shared_ptr<const GenericLinearOperator> A, std::shared_ptr<const GenericLinearOperator> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: const TpetraMatrix& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static std::map<std::string, std::string> methods()
    
        Return a list of available solver methods


    .. cpp:function:: static std::map<std::string, std::string> preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::string parameter_type() const
    
        Return parameter type: "krylov_solver" or "lu_solver"


