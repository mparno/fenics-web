
.. Documentation for the header file dolfin/fem/FiniteElement.h

.. _programmers_reference_cpp_fem_finiteelement:

FiniteElement.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FiniteElement

    This is a wrapper for a UFC finite element (ufc::finite_element).


    .. cpp:function:: FiniteElement(boost::shared_ptr<const ufc::finite_element> element)
    
        Create finite element from UFC finite element (data may be shared)


    .. cpp:function:: std::string signature() const
    
        Return a string identifying the finite element


    .. cpp:function:: ufc::shape cell_shape() const
    
        Return the cell shape


    .. cpp:function:: std::size_t space_dimension() const
    
        Return the dimension of the finite element function space


    .. cpp:function:: std::size_t value_rank() const
    
        Return the rank of the value space


    .. cpp:function:: std::size_t value_dimension(std::size_t i) const
    
        Return the dimension of the value space for axis i


    .. cpp:function:: void evaluate_basis(std::size_t i, double* values, const double* x, const double* vertex_coordinates, int cell_orientation) const
    
        Evaluate basis function i at given point in cell


    .. cpp:function:: void evaluate_basis(std::size_t i, double* values, const double* x, const ufc::cell& c) const
    
        Evaluate basis function i at given point in cell


    .. cpp:function:: void evaluate_basis_all(double* values, const double* x, const double* vertex_coordinates, int cell_orientation) const
    
        Evaluate all basis functions at given point in cell


    .. cpp:function:: void evaluate_basis_all(double* values, const double* x, const double* vertex_coordinates, const ufc::cell& c) const
    
        Evaluate all basis functions at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives(unsigned int i, unsigned int n, double* values, const double* x, const double* vertex_coordinates, int cell_orientation) const
    
        Evaluate order n derivatives of basis function i at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives(unsigned int i, unsigned int n, double* values, const double* x, const ufc::cell& c) const
    
        Evaluate order n derivatives of basis function i at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives_all(unsigned int n, double* values, const double* x, const double* vertex_coordinates, int cell_orientation) const
    
        Evaluate order n derivatives of all basis functions at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives_all(unsigned int n, double* values, const double* x, const ufc::cell& c) const
    
        Evaluate order n derivatives of all basis functions at given point in cell


    .. cpp:function:: double evaluate_dof(std::size_t i, const ufc::function& function, const double* vertex_coordinates, int cell_orientation, const ufc::cell& c) const
    
        Evaluate linear functional for dof i on the function f


    .. cpp:function:: double evaluate_dof(std::size_t i, const ufc::function& function, const ufc::cell& c) const
    
        Evaluate linear functional for dof i on the function f


    .. cpp:function:: void evaluate_dofs(double* values, const ufc::function& f, const double* vertex_coordinates, int cell_orientation, const ufc::cell& c) const
    
        Evaluate linear functionals for all dofs on the function f


    .. cpp:function:: void evaluate_dofs(double* values, const ufc::function& f, const ufc::cell& c) const
    
        Evaluate linear functionals for all dofs on the function f


    .. cpp:function:: void interpolate_vertex_values(double* vertex_values, double* coefficients, int cell_orientation, const ufc::cell& cell) const
    
        Interpolate vertex values from dof values


    .. cpp:function:: void map_from_reference_cell(double* x, const double* xhat, const ufc::cell& c) const
    
        Map coordinate xhat from reference cell to coordinate x in cell


    .. cpp:function:: void map_to_reference_cell(double* xhat, const double* x, const ufc::cell& c) const
    
        Map from coordinate x in cell to coordinate xhat in reference cell


    .. cpp:function:: std::size_t num_sub_elements() const
    
        Return the number of sub elements (for a mixed element)


    .. cpp:function:: std::size_t hash() const
    
        Return simple hash of the signature string


    .. cpp:function:: boost::shared_ptr<const FiniteElement> create_sub_element(std::size_t i) const
    
        Create a new finite element for sub element i (for a mixed element)


    .. cpp:function:: boost::shared_ptr<const FiniteElement> create() const
    
        Create a new class instance


    .. cpp:function:: boost::shared_ptr<const FiniteElement> extract_sub_element(const std::vector<std::size_t>& component) const
    
        Extract sub finite element for component


