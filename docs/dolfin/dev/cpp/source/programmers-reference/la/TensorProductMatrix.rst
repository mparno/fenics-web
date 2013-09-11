
.. Documentation for the header file dolfin/la/TensorProductMatrix.h

.. _programmers_reference_cpp_la_tensorproductmatrix:

TensorProductMatrix.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: TensorProductMatrix

    .. cpp:function:: TensorProductMatrix(std::size_t num_factors)
    
        Create tensor product matrix with given number of factors


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void mult(const TensorProductVector& x, TensorProductVector& y) const
    
        Compute matrix-vector product y = Ax


