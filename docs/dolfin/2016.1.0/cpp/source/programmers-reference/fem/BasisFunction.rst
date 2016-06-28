
.. Documentation for the header file dolfin/fem/BasisFunction.h

.. _programmers_reference_cpp_fem_basisfunction:

BasisFunction.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BasisFunction

    *Parent class(es)*
    
        * :cpp:class:`ufc::function`
        
    This class represents a finite element basis function. It can be
    used for computation of basis function values and derivatives.
    
    Evaluation of basis functions is also possible through the use
    of the functions ``evaluate_basis`` and
    ``evaluate_basis_derivatives`` available in the :cpp:class:`FiniteElement`
    class. The BasisFunction class relies on these functions for
    evaluation but also implements the ufc::function interface which
    allows evaluate_dof to be evaluated for a basis function (on a
    possibly different element).


    .. cpp:function:: BasisFunction(std::size_t index, std::shared_ptr<const FiniteElement> element, const std::vector<double>& coordinate_dofs)
    
        Create basis function with given index on element on given cell
        
        *Arguments*
            index (std::size_t)
                The index of the basis function.
            element (:cpp:class:`FiniteElement`)
                The element to create basis function on.
            cell (ufc::cell)
                The cell.


    .. cpp:function:: void update_index(std::size_t index)
    
        Update the basis function index
        
        *Arguments*
            index (std::size_t)
                The index of the basis function.


    .. cpp:function:: void eval(double* values, const double* x) const
    
        Evaluate basis function at given point
        
        *Arguments*
            values (double)
                The values of the function at the point.
            x (double)
                The coordinates of the point.


    .. cpp:function:: void eval_derivatives(double* values, const double* x, std::size_t n) const
    
        Evaluate all order n derivatives at given point
        
        *Arguments*
            values (double)
                The values of derivatives at the point.
            x (double)
                The coordinates of the point.
            n (std::size_t)
                The order of derivation.


    .. cpp:function:: void evaluate(double* values, const double* coordinates, const ufc::cell& cell) const
    
        Evaluate function at given point in cell
        
        *Arguments*
            values (double)
                The values of the function at the point..
            coordinates (double)
                The coordinates of the point.
            cell (ufc::cell)
                The cell.


