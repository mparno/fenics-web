.. Documentation for the header file dolfin/function/Expression.h

.. _programmers_reference_cpp_function_expression:

Expression.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Expression

    *Parent class*
    
        * :cpp:class:`GenericFunction`
        
    This class represents a user-defined expression. Expressions can
    be used as coefficients in variational forms or interpolated
    into finite element spaces.
    
    An expression is defined by overloading the eval() method. Users
    may choose to overload either a simple version of eval(), in the
    case of expressions only depending on the coordinate x, or an
    optional version for expressions depending on x and mesh data
    like cell indices or facet normals.
    
    The geometric dimension (the size of x) and the value rank and
    dimensions of an expression must supplied as arguments to the
    constructor.

    .. cpp:function:: Expression()
    
        Create scalar expression

    .. cpp:function:: Expression(const Expression& expression)
    
        Copy constructor

    .. cpp:function:: Expression(std::vector<uint> value_shape)
    
        Create tensor-valued expression with given shape

    .. cpp:function:: Expression(uint dim)
    
        Create vector-valued expression with given dimension

    .. cpp:function:: Expression(uint dim0, uint dim1)
    
        Create matrix-valued expression with given dimensions

    .. cpp:function:: uint value_dimension(uint i) const
    
        Return value dimension for given axis

    .. cpp:function:: uint value_rank() const
    
        Return value rank

    .. cpp:function:: void compute_vertex_values(Array<double>& vertex_values,
                                                         const Mesh& mesh) const
    
        Compute values at all mesh vertices

    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x) const
    
        Evaluate expression, must be overloaded by user (simple version)

    .. cpp:function:: void eval(Array<double>& values, const Data& data) const
    
        Evaluate function for given data

    .. cpp:function:: void restrict(double* w,
                                            const FiniteElement& element,
                                            const Cell& dolfin_cell,
                                            const ufc::cell& ufc_cell,
                                            int local_facet) const
    
        Restrict function to local cell (compute expansion coefficients w)

