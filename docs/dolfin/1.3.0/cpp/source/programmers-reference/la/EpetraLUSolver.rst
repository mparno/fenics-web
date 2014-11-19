
.. Documentation for the header file dolfin/la/EpetraLUSolver.h

.. _programmers_reference_cpp_la_epetralusolver:

EpetraLUSolver.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EpetraLUSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLUSolver`
        
    This class implements the direct solution (LU factorization) for
    linear systems of the form Ax = b. It is a wrapper for the LU
    solver of Epetra.


    .. cpp:function:: EpetraLUSolver(std::string method="default")
    
        Constructor


    .. cpp:function:: EpetraLUSolver(boost::shared_ptr<const GenericLinearOperator> A, std::string method="default")
    
        Constructor


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: const GenericLinearOperator& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve(const EpetraMatrix& A, EpetraVector& x, const EpetraVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::size_t solve_transpose(GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve_transpose(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: std::size_t solve_transpose(const EpetraMatrix& A, EpetraVector& x, const EpetraVector& b)
    
        Solve linear system A^Tx = b


    .. cpp:function:: static std::vector<std::pair<std::string, std::string> > methods()
    
        Return a list of available solver methods


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


