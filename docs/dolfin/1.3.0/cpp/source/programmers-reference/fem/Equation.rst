
.. Documentation for the header file dolfin/fem/Equation.h

.. _programmers_reference_cpp_fem_equation:

Equation.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Equation

    This class represents a variational equation lhs == rhs.
    The equation can be either linear or nonlinear:
    
    1. Linear (a == L), in which case a must be a bilinear form
       and L must be a linear form.
    
    2. Nonlinear (F == 0), in which case F must be a linear form.


    .. cpp:function:: Equation(boost::shared_ptr<const Form> a, boost::shared_ptr<const Form> L)
    
        Create equation a == L


    .. cpp:function:: Equation(boost::shared_ptr<const Form> F, int rhs)
    
        Create equation F == 0


    .. cpp:function:: bool is_linear() const
    
        Check whether equation is linear


    .. cpp:function:: boost::shared_ptr<const Form> lhs() const
    
        Return form for left-hand side


    .. cpp:function:: boost::shared_ptr<const Form> rhs() const
    
        Return form for right-hand side


    .. cpp:function:: int rhs_int() const
    
        Return value for right-hand side


