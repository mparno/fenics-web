
.. Documentation for the header file dolfin/la/LinearOperator.h

.. _programmers_reference_cpp_la_linearoperator:

LinearOperator.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LinearOperator

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearOperator`
        
    This class defines an interface for linear operators defined
    only in terms of their action (matrix-vector product) and can be
    used for matrix-free solution of linear systems. The linear
    algebra backend is decided at run-time based on the present
    value of the "linear_algebra_backend" parameter.
    
    To define a linear operator, users need to inherit from this
    class and overload the function mult(x, y) which defines the
    action of the matrix on the vector x as y = Ax.


    .. cpp:function:: LinearOperator()
    
        Create linear operator


    .. cpp:function:: std::size_t size(std::size_t dim) const = 0
    
        Return size of given dimension


    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const = 0
    
        Compute matrix-vector product y = Ax


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: const GenericLinearOperator* instance() const
    
        Return concrete instance / unwrap (const version)


    .. cpp:function:: GenericLinearOperator* instance()
    
        Return concrete instance / unwrap (non-const version)


    .. cpp:function:: boost::shared_ptr<const LinearAlgebraObject> shared_instance() const
    
        Return concrete instance / unwrap (const shared pointer version)


    .. cpp:function:: boost::shared_ptr<LinearAlgebraObject> shared_instance()
    
        Return concrete instance / unwrap (shared pointer version)


