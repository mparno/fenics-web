
.. Documentation for the header file dolfin/function/Constant.h

.. _programmers_reference_cpp_function_constant:

Constant.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Constant

    *Parent class(es)*
    
        * :cpp:class:`Expression`
        
    This class represents a constant-valued expression.


    .. cpp:function:: explicit Constant(double value)
    
        Create scalar constant
        
        *Arguments*
            value (double)
                The scalar to create a Constant object from.
        
        *Example*
            .. code-block:: c++
        
                Constant c(1.0);


    .. cpp:function:: Constant(double value0, double value1)
    
        Create vector constant (dim = 2)
        
        *Arguments*
            value0 (double)
                The first vector element.
            value1 (double)
                The second vector element.
        
        *Example*
            .. code-block:: c++
        
                Constant B(0.0, 1.0);


    .. cpp:function:: Constant(double value0, double value1, double value2)
    
        Create vector constant (dim = 3)
        
        *Arguments*
            value0 (double)
                The first vector element.
            value1 (double)
                The second vector element.
            value2 (double)
                The third vector element.
        
        *Example*
            .. code-block:: c++
        
                Constant T(0.0, 1.0, 0.0);


    .. cpp:function:: explicit Constant(std::vector<double> values)
    
        Create vector-valued constant
        
        *Arguments*
            values (std::vector<double>)
                Values to create a vector-valued constant from.


    .. cpp:function:: Constant(std::vector<std::size_t> value_shape, std::vector<double> values)
    
        Create tensor-valued constant for flattened array of values
        
        *Arguments*
            value_shape (std::vector<std::size_t>)
                Shape of tensor.
            values (std::vector<double>)
                Values to create tensor-valued constant from.


    .. cpp:function:: Constant(const Constant& constant)
    
        Copy constructor
        
        *Arguments*
            constant (:cpp:class:`Constant`)
                Object to be copied.


    .. cpp:function:: const Constant& operator= (const Constant& constant)
    
        Assignment operator
        
        *Arguments*
            constant (:cpp:class:`Constant`)
                Another constant.


    .. cpp:function:: const Constant& operator= (double constant)
    
        Assignment operator
        
        *Arguments*
            constant (double)
                Another constant.


    .. cpp:function:: operator double() const
    
        Cast to double (for scalar constants)
        
        *Returns*
            double
                The scalar value.


    .. cpp:function:: std::vector<double> values() const
    
        Return copy of this Constant's current values
        
        *Returns*
            std::vector<double>
                The vector of scalar values of the constant.


