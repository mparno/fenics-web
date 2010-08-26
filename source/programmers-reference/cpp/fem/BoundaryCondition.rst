.. Documentation for the header file dolfin/fem/BoundaryCondition.h

.. _programmers_reference_cpp_fem_boundarycondition:

BoundaryCondition.h
===================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: BoundaryCondition

    *Parent class*
    
        * :cpp:class:`Variable`
        
        Common base class for boundary conditions

    .. cpp:function:: BoundaryCondition(boost::shared_ptr<const FunctionSpace> V)
    
        Constructor

    .. cpp:function:: BoundaryCondition(const FunctionSpace& V)
    
        Constructor

    .. cpp:function:: boost::shared_ptr<const FunctionSpace> function_space_ptr() const
    
        Return shared pointer to function space

    .. cpp:function:: const FunctionSpace& function_space() const
    
        Return function space

    .. cpp:function:: virtual void apply(GenericMatrix& A) const = 0
    
        Apply boundary condition to a matrix

    .. cpp:function:: virtual void apply(GenericMatrix& A, GenericVector& b) const = 0
    
        Apply boundary condition to a linear system

    .. cpp:function:: virtual void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const = 0
    
        Apply boundary condition to a linear system for a nonlinear problem

    .. cpp:function:: virtual void apply(GenericVector& b) const = 0
    
        Apply boundary condition to a vector

    .. cpp:function:: virtual void apply(GenericVector& b, const GenericVector& x) const = 0
    
        Apply boundary condition to a vector for a nonlinear problem

    .. cpp:function:: virtual ~BoundaryCondition()
    
        Destructor

.. cpp:class:: LocalData

