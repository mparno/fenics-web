
.. Documentation for the header file dolfin/function/FunctionAXPY.h

.. _programmers_reference_cpp_function_functionaxpy:

FunctionAXPY.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FunctionAXPY

    This class represents a linear combination of functions. It is
    mostly used as an intermediate class for operations such as u =
    3*u0 + 4*u1; where the rhs generates an FunctionAXPY.


    .. cpp:function:: FunctionAXPY(const Function& func, double scalar)
    
        Constructor


    .. cpp:function:: FunctionAXPY(const FunctionAXPY& axpy, double scalar)
    
        Constructor


    .. cpp:function:: FunctionAXPY(const Function& func0, const Function& func1, Direction direction)
    
        Constructor


    .. cpp:function:: FunctionAXPY(const FunctionAXPY& axpy, const Function& func, Direction direction)
    
        Constructor


    .. cpp:function:: FunctionAXPY(const FunctionAXPY& axpy0, const FunctionAXPY& axpy1, Direction direction)
    
        Constructor


    .. cpp:function:: FunctionAXPY(std::vector<std::pair<double, const Function*> > pairs)
    
        Constructor


    .. cpp:function:: FunctionAXPY(const FunctionAXPY& axpy)
    
        Copy constructor


    .. cpp:function:: FunctionAXPY operator+(const Function& func) const
    
        Addition operator


    .. cpp:function:: FunctionAXPY operator+(const FunctionAXPY& axpy) const
    
        Addition operator


    .. cpp:function:: FunctionAXPY operator-(const Function& func) const
    
        Substraction operator


    .. cpp:function:: FunctionAXPY operator-(const FunctionAXPY& axpy) const
    
        Substraction operator


    .. cpp:function:: FunctionAXPY operator*(double scale) const
    
        Scale operator


    .. cpp:function:: FunctionAXPY operator/(double scale) const
    
        Scale operator


    .. cpp:function:: const std::vector<std::pair<double, const Function*> >& pairs() const
    
        Return the scalar and Function pairs


    .. cpp:function:: void _register(const FunctionAXPY& axpy0, double scale)
    
        Register another AXPY object


