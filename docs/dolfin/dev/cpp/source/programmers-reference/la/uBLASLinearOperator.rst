
.. Documentation for the header file dolfin/la/uBLASLinearOperator.h

.. _programmers_reference_cpp_la_ublaslinearoperator:

uBLASLinearOperator.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASLinearOperator

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearOperator`
        
    .. cpp:function:: uBLASLinearOperator()
    
        Constructor


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const
    
        Compute matrix-vector product y = Ax


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


