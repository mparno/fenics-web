
.. Documentation for the header file dolfin/la/GenericLinearOperator.h

.. _programmers_reference_cpp_la_genericlinearoperator:

GenericLinearOperator.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericLinearOperator

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraObject`
        
    This class defines a common interface for linear operators,
    including actual matrices (class :cpp:class:`GenericMatrix`) and linear
    operators only defined in terms of their action on vectors.
    
    This class is used internally by DOLFIN to define a class
    hierarchy of backend independent linear operators and solvers.
    Users should not interface to this class directly but instead
    use the :cpp:class:`LinearOperator` class.


    .. cpp:function:: std::size_t size(std::size_t dim) const = 0
    
        Return size of given dimension


    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const = 0
    
        Compute matrix-vector product y = Ax


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


