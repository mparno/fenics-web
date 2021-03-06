
.. Documentation for the header file dolfin/la/ITLKrylovSolver.h

.. _programmers_reference_cpp_la_itlkrylovsolver:

ITLKrylovSolver.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: ITLKrylovSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearSolver`
        
    This class implements Krylov methods for linear systems
    of the form Ax = b. It is a wrapper for the Krylov solvers
    of ITL.


    .. cpp:function:: ITLKrylovSolver(std::string method = "default", std::string preconditioner = "default")
    
        Create Krylov solver for a particular method and preconditioner


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operators(const boost::shared_ptr<const GenericMatrix> A, const boost::shared_ptr<const GenericMatrix> P)
    
        Set operator (matrix) and preconditioner matrix


    .. cpp:function:: const GenericMatrix& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(MTL4Vector& x, const MTL4Vector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b and return number of iterations


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available methods


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > preconditioners()
    
        Return a list available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


