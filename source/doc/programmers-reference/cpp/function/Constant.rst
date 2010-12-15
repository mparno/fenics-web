.. Documentation for the header file dolfin/function/Constant.h

.. _programmers_reference_cpp_function_constant:

Constant.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Constant

    *Parent class*
    
        * :cpp:class:`Expression`
        
    This class represents a constant-valued expression.

    .. cpp:function:: Constant(double value)
    
        Create scalar constant

    .. cpp:function:: Constant(double value0, double value1)
    
        Create vector constant (dim = 2)

    .. cpp:function:: Constant(double value0, double value1, double value2)
    
        Create vector constant (dim = 3)

    .. cpp:function:: Constant(std::vector<double> values)
    
        Create vector-valued constant

    .. cpp:function:: Constant(std::vector<uint> value_shape, std::vector<double> values)
    
        Create tensor-valued constant for flattened array of values

    .. cpp:function:: Constant(const Constant& constant)
    
        Copy constructor

    .. cpp:function:: const Constant& operator= (const Constant& constant)
    
        Assignment operator

    .. cpp:function:: const Constant& operator= (double constant)
    
        Assignment operator

    .. cpp:function:: operator double() const
    
        Cast to double (for scalar constants)

