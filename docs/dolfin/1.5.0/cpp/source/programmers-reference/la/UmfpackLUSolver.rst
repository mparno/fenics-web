
.. Documentation for the header file dolfin/la/UmfpackLUSolver.h

.. _programmers_reference_cpp_la_umfpacklusolver:

UmfpackLUSolver.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UmfpackLUSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLUSolver`
        
    This class implements the direct solution (LU factorization) of
    linear systems of the form Ax = b using UMFPACK
    (http://www.cise.ufl.edu/research/sparse/umfpack/) if installed.


    .. cpp:function:: UmfpackLUSolver()
    
        Constructor


    .. cpp:function:: UmfpackLUSolver(std::shared_ptr<const GenericLinearOperator> A)
    
        Constructor


    .. cpp:function:: void set_operator(std::shared_ptr<const GenericLinearOperator> A)
    
        Set operator (matrix)


    .. cpp:function:: const GenericLinearOperator& get_operator() const
    
        Return the operator (matrix)


    .. cpp:function:: std::size_t solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b for a sparse matrix using UMFPACK
        if installed


    .. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


    .. cpp:function:: static void umfpack_check_status(long int status, std::string function)
    
        Check status flag returned by an UMFPACK function


