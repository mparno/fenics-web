
.. Documentation for the header file dolfin/la/PETScLinearOperator.h

.. _programmers_reference_cpp_la_petsclinearoperator:

PETScLinearOperator.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScLinearOperator

    *Parent class(es)*
    
        * :cpp:class:`PETScBaseMatrix`
        
    .. cpp:function:: PETScLinearOperator()
    
        Constructor


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Return size of given dimension


    .. cpp:function:: void mult(const GenericVector& x, GenericVector& y) const
    
        Compute matrix-vector product y = Ax


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: const GenericLinearOperator* wrapper() const
    
        Return pointer to wrapper (const version)


    .. cpp:function:: GenericLinearOperator* wrapper()
    
        Return pointer to wrapper (const version)


