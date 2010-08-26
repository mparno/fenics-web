.. Documentation for the header file dolfin/function/Function.h

.. _programmers_reference_cpp_function_Mesh:

Function.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Function

    *Parent class*
    
        * :cpp:class:`GenericFunction`
        
        This class represents a function u_h in a finite element
        function space V_h, given by
        
        u_h = sum_i U_i phi_i
        
        where {phi_i}_i is a basis for V_h, and U is a vector of
        expansion coefficients for u_h.

    .. cpp:function:: Function& operator[] (uint i) const
    
        Extract sub-function

    .. cpp:function:: Function(boost::shared_ptr<const FunctionSpace> V,
                               GenericVector& x)
    
        Create function on given function space with a given vector (used by
        Python interface)

    .. cpp:function:: Function(boost::shared_ptr<const FunctionSpace> V,
                               boost::shared_ptr<GenericVector> x)
    
        Create function on given function space with a given vector
        (shared data)

    .. cpp:function:: Function(boost::shared_ptr<const FunctionSpace> V,
                               std::string filename)
    
        Create function from vector of dofs stored to file (shared data)

    .. cpp:function:: Function(const Function& v)
    
        Copy constructor

    .. cpp:function:: Function(const Function& v, uint i)
    
        Sub-function constructor with shallow copy of vector (used in Python
        interface)

    .. cpp:function:: Function(const FunctionSpace& V,
                               GenericVector& x)
    
        Create function on given function space with a given vector

    .. cpp:function:: Function(const FunctionSpace& V,
                               std::string filename)
    
        Create function from vector of dofs stored to file

    .. cpp:function:: GenericVector& vector()
    
        Return vector of expansion coefficients (non-const version)

    .. cpp:function:: bool in(const FunctionSpace& V) const
    
        Check if function is a member of the given function space

    .. cpp:function:: boost::shared_ptr<const FunctionSpace> function_space_ptr() const
    
        Return shared pointer to function space

    .. cpp:function:: const Function& operator= (const Expression& v)
    
        Assignment from expression using interpolation

    .. cpp:function:: const Function& operator= (const Function& v)
    
        Assignment from function

    .. cpp:function:: const FunctionSpace& function_space() const
    
        Return function space

    .. cpp:function:: const GenericVector& vector() const
    
        Return vector of expansion coefficients (const version)

    .. cpp:function:: explicit Function(boost::shared_ptr<const FunctionSpace> V)
    
        Create function on given function space (shared data)

    .. cpp:function:: explicit Function(const FunctionSpace& V)
    
        Create function on given function space

    .. cpp:function:: uint geometric_dimension() const
    
        Return geometric dimension

    .. cpp:function:: virtual uint value_dimension(uint i) const
    
        Return value dimension for given axis

    .. cpp:function:: virtual uint value_rank() const
    
        Return value rank

    .. cpp:function:: virtual void compute_vertex_values(Array<double>& vertex_values,
                                                         const Mesh& mesh) const
    
        Compute values at all mesh vertices

    .. cpp:function:: virtual void eval(Array<double>& values, const Data& data) const
    
        Evaluate function for given data

    .. cpp:function:: virtual void gather() const
    
        Collect off-process coefficients to prepare for interpolation

    .. cpp:function:: virtual void restrict(double* w,
                                            const FiniteElement& element,
                                            const Cell& dolfin_cell,
                                            const ufc::cell& ufc_cell,
                                            int local_facet) const
    
        Restrict function to local cell (compute expansion coefficients w)

    .. cpp:function:: virtual ~Function()
    
        Destructor

    .. cpp:function:: void eval(Array<double>& values,
                                const Array<double>& x,
                                const Cell& dolfin_cell,
                                const ufc::cell& ufc_cell) const
    
        Evaluate function for given coordinate in given cell

    .. cpp:function:: void eval(Array<double>& values, const Array<double>& x) const
    
        Evaluate function for given coordinate

    .. cpp:function:: void extrapolate(const Function& v)
    
        Extrapolate function (from a possibly lower-degree function space)

    .. cpp:function:: void interpolate(const GenericFunction& v)
    
        Interpolate function (possibly non-matching meshes)

.. cpp:class:: LocalScratch

.. cpp:class:: GatherScratch

