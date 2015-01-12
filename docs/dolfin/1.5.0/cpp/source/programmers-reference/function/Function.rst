
.. Documentation for the header file dolfin/function/Function.h

.. _programmers_reference_cpp_function_function:

Function.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Function

    *Parent class(es)*
    
        * :cpp:class:`GenericFunction`
        
        * :cpp:class:`Hierarchical<Function>`
        
    This class represents a function :math:`u_h` in a finite
    element function space :math:`V_h`, given by
    
    .. math::
    
        u_h = \sum_{i=1}^{n} U_i \phi_i
    
    where :math:`\{\phi_i\}_{i=1}^{n}` is a basis for :math:`V_h`,
    and :math:`U` is a vector of expansion coefficients for :math:`u_h`.


    .. cpp:function:: explicit Function(const FunctionSpace& V)
    
        Create function on given function space
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
        
        *Example*
            .. code-block:: c++
        
                Function u(V);
        


    .. cpp:function:: explicit Function(std::shared_ptr<const FunctionSpace> V)
    
        Create function on given function space (shared data)
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: Function(std::shared_ptr<const FunctionSpace> V, std::shared_ptr<GenericVector> x)
    
        Create function on given function space with a given vector
        (shared data)
        
        *Warning: This constructor is intended for internal library use only*
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            x (:cpp:class:`GenericVector`)
                The vector.


    .. cpp:function:: Function(const FunctionSpace& V, std::string filename)
    
        Create function from vector of dofs stored to file
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            filename_vector (std::string)
                The name of the file containing the vector.
            filename_dofdata (std::string)
                The name of the file containing the dofmap data.


    .. cpp:function:: Function(std::shared_ptr<const FunctionSpace> V, std::string filename)
    
        Create function from vector of dofs stored to file (shared data)
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            filename_dofdata (std::string)
                The name of the file containing the dofmap data.


    .. cpp:function:: Function(const Function& v)
    
        Copy constructor
        
        *Arguments*
            v (:cpp:class:`Function`)
                The object to be copied.


    .. cpp:function:: Function(const Function& v, std::size_t i)
    
        Sub-function constructor with shallow copy of vector (used in Python
        interface)
        
        *Arguments*
            v (:cpp:class:`Function`)
                The function to be copied.
            i (std::size_t)
                Index of subfunction.
        


    .. cpp:function:: const Function& operator= (const Function& v)
    
        Assignment from function
        
        *Arguments*
            v (:cpp:class:`Function`)
                Another function.


    .. cpp:function:: const Function& operator= (const Expression& v)
    
        Assignment from expression using interpolation
        
        *Arguments*
            v (:cpp:class:`Expression`)
                The expression.


    .. cpp:function:: void operator=(const FunctionAXPY& axpy)
    
        Assignment from linear combination of function
        
        *Arguments*
            v (:cpp:class:`FunctionAXPY`)
                A linear combination of other Functions


    .. cpp:function:: Function& operator[] (std::size_t i) const
    
        Extract subfunction
        
        *Arguments*
            i (std::size_t)
                Index of subfunction.
        *Returns*
            :cpp:class:`Function`
                The subfunction.


    .. cpp:function:: FunctionAXPY operator+(const Function& other) const
    
        Add operator with other function
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: FunctionAXPY operator+(const FunctionAXPY& axpy) const
    
        Add operator with other linear combination of functions
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: FunctionAXPY operator-(const Function& other) const
    
        Subtraction operator with other function
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: FunctionAXPY operator-(const FunctionAXPY& axpy) const
    
        Subtraction operator with other linear combination of functions
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: FunctionAXPY operator*(double scalar) const
    
        Scale operator
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: FunctionAXPY operator/(double scalar) const
    
        Scale operator
        
        *Returns*
            :cpp:class:`FunctionAXPY`
                Return a linear combination of Functions


    .. cpp:function:: std::shared_ptr<const FunctionSpace> function_space() const
    
        Return shared pointer to function space
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Return the shared pointer.


    .. cpp:function:: std::shared_ptr<GenericVector> vector()
    
        Return vector of expansion coefficients (non-const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients.


    .. cpp:function:: std::shared_ptr<const GenericVector> vector() const
    
        Return vector of expansion coefficients (const version)
        
        *Returns*
            :cpp:class:`GenericVector`
                The vector of expansion coefficients (const).


    .. cpp:function:: bool in(const FunctionSpace& V) const
    
        Check if function is a member of the given function space
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
        
        *Returns*
            bool
                True if the function is in the function space.


    .. cpp:function:: std::size_t geometric_dimension() const
    
        Return geometric dimension
        
        *Returns*
            std::size_t
                The geometric dimension.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x) const
    
        Evaluate function at given coordinates
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values.
            x (:cpp:class:`Array` <double>)
                The coordinates.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const Cell& dolfin_cell, const ufc::cell& ufc_cell) const
    
        Evaluate function at given coordinates in given cell
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values.
            x (:cpp:class:`Array` <double>)
                The coordinates.
            dolfin_cell (:cpp:class:`Cell`)
                The cell.
            ufc_cell (ufc::cell)
                The ufc::cell.


    .. cpp:function:: void interpolate(const GenericFunction& v)
    
        Interpolate function (on possibly non-matching meshes)
        
        *Arguments*
            v (:cpp:class:`GenericFunction`)
                The function to be interpolated.


    .. cpp:function:: void extrapolate(const Function& v)
    
        Extrapolate function (from a possibly lower-degree function space)
        
        *Arguments*
            v (:cpp:class:`Function`)
                The function to be extrapolated.


    .. cpp:function:: std::size_t value_rank() const
    
        Return value rank
        
        *Returns*
            std::size_t
                The value rank.


    .. cpp:function:: std::size_t value_dimension(std::size_t i) const
    
        Return value dimension for given axis
        
        *Arguments*
            i (std::size_t)
                The index of the axis.
        
        *Returns*
            std::size_t
                The value dimension.


    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
    
        Evaluate at given point in given cell
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values at the point.
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.
            cell (ufc::cell)
                The cell which contains the given point.


    .. cpp:function:: void non_matching_eval(Array<double>& values, const Array<double>& x, const ufc::cell& ufc_cell) const
    
        Evaluate function for given data (non-matching meshes)
        
        *Arguments*
            values (:cpp:class:`Array` <double>)
                The values at the point.
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.
            cell (ufc::cell)
                The cell.


    .. cpp:function:: void restrict(double* w, const FiniteElement& element, const Cell& dolfin_cell, const double* vertex_coordinates, const ufc::cell& ufc_cell) const
    
        Restrict function to local cell (compute expansion coefficients w)
        
        *Arguments*
            w (list of doubles)
                Expansion coefficients.
            element (:cpp:class:`FiniteElement`)
                The element.
            dolfin_cell (:cpp:class:`Cell`)
                The cell.
            ufc_cell (ufc::cell).
                The ufc::cell.


    .. cpp:function:: void compute_vertex_values(std::vector<double>& vertex_values, const Mesh& mesh) const
    
        Compute values at all mesh vertices
        
        *Arguments*
            vertex_values (:cpp:class:`Array` <double>)
                The values at all vertices.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: void compute_vertex_values(std::vector<double>& vertex_values)
    
        Compute values at all mesh vertices
        
        *Arguments*
            vertex_values (:cpp:class:`Array` <double>)
                The values at all vertices.


