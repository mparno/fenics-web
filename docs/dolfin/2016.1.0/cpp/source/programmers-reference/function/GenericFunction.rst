
.. Documentation for the header file dolfin/function/GenericFunction.h

.. _programmers_reference_cpp_function_genericfunction:

GenericFunction.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericFunction

    *Parent class(es)*
    
        * :cpp:class:`ufc::function`
        
        * :cpp:class:`Variable`
        
    This is a common base class for functions. Functions can be
    evaluated at a given point and they can be restricted to a given
    cell in a finite element mesh. This functionality is implemented
    by sub-classes that implement the eval() and restrict() functions.
    
    DOLFIN provides two implementations of the GenericFunction
    interface in the form of the classes Function and Expression.
    
    Sub-classes may optionally implement the update() function that
    will be called prior to restriction when running in parallel.


    .. cpp:function:: GenericFunction()
    
        Constructor


    .. cpp:function:: std::size_t value_rank() const = 0
    
        Return value rank


    .. cpp:function:: std::size_t value_dimension(std::size_t i) const = 0
    
        Return value dimension for given axis


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Evaluate at given point in given cell


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x) const
    
        Evaluate at given point


    .. cpp:function:: void restrict(double* w, const FiniteElement& element, const Cell& dolfin_cell, const double* coordinate_dofs, const ufc::cell& ufc_cell) const = 0
    
        Restrict function to local cell (compute expansion coefficients w)


    .. cpp:function:: void compute_vertex_values(std::vector<double>& vertex_values, const Mesh& mesh) const = 0
    
        Compute values at all mesh vertices


    .. cpp:function:: void update() const
    
        Update off-process ghost coefficients


    .. cpp:function:: double operator() (double x) const
    
        Evaluation at given point (scalar function)


    .. cpp:function:: double operator() (double x, double y) const
    
        Evaluation at given point (scalar function)


    .. cpp:function:: double operator() (double x, double y, double z) const
    
        Evaluation at given point (scalar function)


    .. cpp:function:: double operator() (const Point& p) const
    
        Evaluation at given point (scalar function)


    .. cpp:function:: void operator() (Array<double>& values, double x) const
    
        Evaluation at given point (vector-valued function)


    .. cpp:function:: void operator() (Array<double>& values, double x, double y) const
    
        Evaluation at given point (vector-valued function)


    .. cpp:function:: void operator() (Array<double>& values, double x, double y, double z) const
    
        Evaluation at given point (vector-valued function)


    .. cpp:function:: void operator() (Array<double>& values, const Point& p) const
    
        Evaluation at given point (vector-valued function)


    .. cpp:function:: std::size_t value_size() const
    
        Evaluation at given point
        Return value size (product of value dimensions)


    .. cpp:function:: void evaluate(double* values, const double* coordinates, const ufc::cell& cell) const
    
        Evaluate function at given point in cell


