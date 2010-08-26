.. Documentation for the header file dolfin/la/solve.h

.. _programmers_reference_cpp_la_Mesh:

solve.h
=======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

    .. cpp:function:: double normalize(GenericVector& x, std::string normalization_type = "average")
    
        Normalize vector according to given normalization type

    .. cpp:function:: double residual(const GenericMatrix& A, const GenericVector& x, const GenericVector& b)
    
        Compute residual ||Ax - b||

    .. cpp:function:: void solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b,
                                 std::string solver_type = "lu", std::string pc_type = "default")
    
        Solve linear system Ax = b

