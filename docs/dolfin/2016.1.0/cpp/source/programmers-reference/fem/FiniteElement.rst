
.. Documentation for the header file dolfin/fem/FiniteElement.h

.. _programmers_reference_cpp_fem_finiteelement:

FiniteElement.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FiniteElement

    This is a wrapper for a UFC finite element (ufc::finite_element).


    .. cpp:function:: FiniteElement(std::shared_ptr<const ufc::finite_element> element)
    
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


    .. cpp:function:: void evaluate_basis(std::size_t i, double* values, const double* x, const double* coordinate_dofs, int cell_orientation) const
    
        Evaluate basis function i at given point in cell


    .. cpp:function:: void evaluate_basis_all(double* values, const double* x, const double* coordinate_dofs, int cell_orientation) const
    
        Evaluate all basis functions at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives(unsigned int i, unsigned int n, double* values, const double* x, const double* coordinate_dofs, int cell_orientation) const
    
        Evaluate order n derivatives of basis function i at given point in cell


    .. cpp:function:: void evaluate_basis_derivatives_all(unsigned int n, double* values, const double* x, const double* coordinate_dofs, int cell_orientation) const
    
        Evaluate order n derivatives of all basis functions at given
        point in cell


    .. cpp:function:: double evaluate_dof(std::size_t i, const ufc::function& function, const double* coordinate_dofs, int cell_orientation, const ufc::cell& c) const
    
        Evaluate linear functional for dof i on the function f


    .. cpp:function:: void evaluate_dofs(double* values, const ufc::function& f, const double* coordinate_dofs, int cell_orientation, const ufc::cell& c) const
    
        Evaluate linear functionals for all dofs on the function f


    .. cpp:function:: void interpolate_vertex_values(double* vertex_values, double* coefficients, const double* coordinate_dofs, int cell_orientation, const ufc::cell& cell) const
    
        Interpolate vertex values from dof values


    .. cpp:function:: void tabulate_dof_coordinates(boost::multi_array<double, 2>& coordinates, const std::vector<double>& coordinate_dofs, const Cell& cell) const
    
        Tabulate the coordinates of all dofs on an element
        
        *Arguments*
            coordinates (boost::multi_array<double, 2>)
                The coordinates of all dofs on a cell.
            coordinate_dofs (std::vector<double>)
                The cell coordinates
            cell (Cell)
                The cell.


    .. cpp:function:: std::size_t num_sub_elements() const
    
        Return the number of sub elements (for a mixed element)


    .. cpp:function:: std::size_t hash() const
    
        Return simple hash of the signature string


    .. cpp:function:: std::shared_ptr<const FiniteElement> create_sub_element(std::size_t i) const
    
        Create a new finite element for sub element i (for a mixed element)


    .. cpp:function:: std::shared_ptr<const FiniteElement> create() const
    
        Create a new class instance


    .. cpp:function:: std::shared_ptr<const FiniteElement> extract_sub_element(const std::vector<std::size_t>& component) const
    
        Extract sub finite element for component


    .. cpp:function:: std::shared_ptr<const ufc::finite_element> ufc_element() const
    
        Return underlying UFC element. Intended for libray usage only
        and may change.


