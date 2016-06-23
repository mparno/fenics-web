
.. Documentation for the header file dolfin/function/Expression.h

.. _programmers_reference_cpp_function_expression:

Expression.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Expression

    *Parent class(es)*
    
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
    
        Create scalar expression.


    .. cpp:function:: explicit Expression(std::size_t dim)
    
        Create vector-valued expression with given dimension.
        
        *Arguments*
            dim (std::size_t)
                Dimension of the vector-valued expression.


    .. cpp:function:: Expression(std::size_t dim0, std::size_t dim1)
    
        Create matrix-valued expression with given dimensions.
        
        *Arguments*
            dim0 (std::size_t)
                Dimension (rows).
            dim1 (std::size_t)
                Dimension (columns).


    .. cpp:function:: explicit Expression(std::vector<std::size_t> value_shape)
    
        Create tensor-valued expression with given shape.
        
        *Arguments*
            value_shape (std::vector<std::size_t>)
                Shape of expression.


    .. cpp:function:: Expression(const Expression& expression)
    
        Copy constructor
        
        *Arguments*
            expression (:cpp:class:`Expression`)
                Object to be copied.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Note: The reimplementation of eval is needed for the Python interface.
        Evaluate at given point in given cell.
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values at the point.
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.
            cell (ufc::cell)
                The cell which contains the given point.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x) const
    
        Evaluate at given point.
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values at the point.
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.


    .. cpp:function:: std::size_t value_rank() const
    
        Return value rank.
        
        *Returns*
            std::size_t
                The value rank.


    .. cpp:function:: std::size_t value_dimension(std::size_t i) const
    
        Return value dimension for given axis.
        
        *Arguments*
            i (std::size_t)
                Integer denoting the axis to use.
        
        *Returns*
            std::size_t
                The value dimension (for the given axis).


    .. cpp:function:: void restrict(double* w, const FiniteElement& element, const Cell& dolfin_cell, const double* coordinate_dofs, const ufc::cell& ufc_cell) const
    
        Restrict function to local cell (compute expansion coefficients w).
        
        *Arguments*
            w (list of doubles)
                Expansion coefficients.
            element (:cpp:class:`FiniteElement`)
                The element.
            dolfin_cell (:cpp:class:`Cell`)
                The cell.
            ufc_cell (ufc::cell)
                The ufc::cell.


    .. cpp:function:: void compute_vertex_values(std::vector<double>& vertex_values, const Mesh& mesh) const
    
        Compute values at all mesh vertices.
        
        *Arguments*
            vertex_values (:cpp:class:`Array` <double>)
                The values at all vertices.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: std::shared_ptr<const FunctionSpace> function_space() const
    
        Return shared pointer to function space (NULL)
        Expression does not have a FunctionSpace
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Return the shared pointer.


