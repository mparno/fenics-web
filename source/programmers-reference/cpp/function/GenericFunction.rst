.. Documentation for the header file dolfin/function/GenericFunction.h

.. _programmers_reference_cpp_function_Mesh:

GenericFunction.h
=================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: GenericFunction

    *Parent class*
    
        * :cpp:class:`ufc::function,`
        
        This is a common base class for functions. Functions can be
        evaluated at a given point and they can be restricted to a given
        cell in a finite element mesh. This functionality is implemented
        by sub-classes that implement the eval() and restrict() functions.
        
        DOLFIN provides two implementations of the GenericFunction
        interface in the form of the classes Function and Expression.
        
        Sub-classes may optionally implement the gather() function that
        will be called prior to restriction when running in parallel.

    .. cpp:function:: GenericFunction()
    
        Constructor

    .. cpp:function:: uint value_size() const
    
        Return value size (product of value dimensions)

    .. cpp:function:: virtual uint value_dimension(uint i) const = 0
    
        Return value dimension for given axis

    .. cpp:function:: virtual uint value_rank() const = 0
    
        Return value rank

    .. cpp:function:: virtual void compute_vertex_values(Array<double>& vertex_values,
                                                         const Mesh& mesh) const = 0
    
        Compute values at all mesh vertices

    .. cpp:function:: virtual void eval(Array<double>& values, const Data& data) const = 0
    
        Evaluate function for given data

    .. cpp:function:: virtual void evaluate(double* values,
                                            const double* coordinates,
                                            const ufc::cell& cell) const
    
        Evaluate function at given point in cell

    .. cpp:function:: virtual void gather() const
    
        Collect off-process coefficients to prepare for interpolation

    .. cpp:function:: virtual void restrict(double* w,
                                            const FiniteElement& element,
                                            const Cell& dolfin_cell,
                                            const ufc::cell& ufc_cell,
                                            int local_facet) const = 0
    
        Restrict function to local cell (compute expansion coefficients w)

    .. cpp:function:: virtual ~GenericFunction()
    
        Destructor

    .. cpp:function:: void restrict(double* w,
                                    const FiniteElement& element,
                                    const Cell& dolfin_cell,
                                    const ufc::cell& ufc_cell) const
    
        Convenience function for restriction when facet is unknown

