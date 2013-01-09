
.. Documentation for the header file dolfin/fem/BoundaryCondition.h

.. _programmers_reference_cpp_fem_boundarycondition:

BoundaryCondition.h
===================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BoundaryCondition

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    Common base class for boundary conditions


    .. cpp:function:: BoundaryCondition(const FunctionSpace& V)
    
        Constructor


    .. cpp:function:: BoundaryCondition(boost::shared_ptr<const FunctionSpace> V)
    
        Constructor


    .. cpp:function:: void apply(GenericMatrix& A) const = 0
    
        Apply boundary condition to a matrix


    .. cpp:function:: void apply(GenericVector& b) const = 0
    
        Apply boundary condition to a vector


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b) const = 0
    
        Apply boundary condition to a linear system


    .. cpp:function:: void apply(GenericVector& b, const GenericVector& x) const = 0
    
        Apply boundary condition to a vector for a nonlinear problem


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const = 0
    
        Apply boundary condition to a linear system for a nonlinear problem


    .. cpp:function:: boost::shared_ptr<const FunctionSpace> function_space() const
    
        Return shared pointer to function space


.. cpp:class:: LocalData

