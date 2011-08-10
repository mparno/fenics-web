
.. Documentation for the header file dolfin/la/PETScLUSolver.h

.. _programmers_reference_cpp_la_petsclusolver:

PETScLUSolver.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScLUSolver

    *Parent class(es)*
    
        * :cpp:class:`GenericLUSolver`
        
        * :cpp:class:`PETScObject`
        
    This class implements the direct solution (LU factorization) for
    linear systems of the form Ax = b. It is a wrapper for the LU
    solver of PETSc.


    .. cpp:function:: PETScLUSolver(std::string lu_package="default")
    
        Constructor


    .. cpp:function:: PETScLUSolver(boost::shared_ptr<const PETScMatrix> A, std::string lu_package="default")
    
        Constructor


    .. cpp:function:: void set_operator(const boost::shared_ptr<const GenericMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: void set_operator(const boost::shared_ptr<const PETScMatrix> A)
    
        Set operator (matrix)


    .. cpp:function:: const GenericMatrix& get_operator() const
    
        Get operator (matrix)


    .. cpp:function:: uint solve(GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: uint solve(const PETScMatrix& A, PETScVector& x, const PETScVector& b)
    
        Solve linear system Ax = b


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: boost::shared_ptr<KSP> ksp() const
    
        Return PETSc KSP pointer


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


