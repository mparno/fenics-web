
.. Documentation for the header file dolfin/la/EigenLUSolver.h

.. _programmers_reference_cpp_la_eigenlusolver:

EigenLUSolver.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EigenLUSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLUSolver`
        
    This class implements the direct solution (LU factorization) for
    linear systems of the form Ax = b.


    .. cpp:function:: EigenLUSolver(std::string method="default")
    
        Constructor


    .. cpp:function:: EigenLUSolver(std::shared_ptr<const EigenMatrix> A, std::string method="default")
    
        Constructor


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operator(std::shared_ptr<const EigenMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: const GenericLinearOperator& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b, bool transpose)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(const EigenMatrix& A, EigenVector& x, const EigenVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve_transpose(GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve_transpose(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve_transpose(const EigenMatrix& A, EigenVector& x, const EigenVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: static std::map<std::string, std::string> methods()
    
        Return a list of available solver methods


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


