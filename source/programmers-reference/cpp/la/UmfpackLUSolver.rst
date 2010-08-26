.. Documentation for the header file dolfin/la/UmfpackLUSolver.h

.. _programmers_reference_cpp_la_Mesh:

UmfpackLUSolver.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: class GenericVector
    
        Forward declarations

.. cpp:class:: UmfpackLUSolver

    *Parent class*
    
        * :cpp:class:`GenericLUSolver`
        
        This class implements the direct solution (LU factorization) of
        linear systems of the form Ax = b using UMFPACK
        (http://www.cise.ufl.edu/research/sparse/umfpack/) if installed.

    .. cpp:function:: UmfpackLUSolver()
    
        Constructor

    .. cpp:function:: UmfpackLUSolver(boost::shared_ptr<const GenericMatrix> A)
    
        Constructor

    .. cpp:function:: UmfpackLUSolver(const GenericMatrix& A)
    
        Constructor

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: static void umfpack_check_status(long int status, std::string function)
    
        Check status flag returned by an UMFPACK function

    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b for a sparse matrix using UMFPACK if installed

    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system

    .. cpp:function:: uint solve_factorized(GenericVector& x, const GenericVector& b) const
    
        Solve factorized system (UMFPACK).

    .. cpp:function:: void numeric_factorize()
    
        LU factorisation

    .. cpp:function:: void set_operator(const GenericMatrix& A)
    
        Set operator (matrix)

    .. cpp:function:: ~UmfpackLUSolver()
    
        Destructor

