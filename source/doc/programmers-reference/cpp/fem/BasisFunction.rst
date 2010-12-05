.. Documentation for the header file dolfin/fem/BasisFunction.h

.. _programmers_reference_cpp_fem_basisfunction:

BasisFunction.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: BasisFunction

    *Parent class*
    
        * :cpp:class:`ufc::function`
        
    This class represents a finite element basis function. It can be
    used for computation of basis function values and derivatives.
    
    Evaluation of basis functions is also possible through the use
    of the functions evaluate_basis and evaluate_basis_derivatives
    available in the FiniteElement class. The BasisFunction class
    relies on these functions for evaluation but also implements the
    ufc::function interface which allows evaluate_dof to be
    evaluated for a basis function (on a possibly different
    element).

    .. cpp:function:: BasisFunction(uint index, const FiniteElement& element, const ufc::cell& cell)
    
        Create basis function with given index on element on given cell

    .. cpp:function:: void eval(double* values, const double* x) const
    
        Evaluate basis function at given point

    .. cpp:function:: void eval_derivatives(double* values, const double* x, uint n) const
    
        Evaluate all order n derivatives at given point

    .. cpp:function:: void evaluate(double* values, const double* coordinates, const ufc::cell& cell) const
    
        Evaluate function at given point in cell

