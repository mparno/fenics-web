.. Documentation for the header file dolfin/fem/FiniteElement.h

.. _programmers_reference_cpp_fem_finiteelement:

FiniteElement.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: FiniteElement

    This is a wrapper for a UFC finite element (ufc::finite_element).

    .. cpp:function:: FiniteElement(boost::shared_ptr<const ufc::finite_element> element)
    
        Create finite element from UFC finite element (data may be shared)

    .. cpp:function:: uint hash() const
    
        Return simple hash of the signature string

    .. cpp:function:: boost::shared_ptr<const FiniteElement> create_sub_element(uint i) const
    
        Create sub element

    .. cpp:function:: boost::shared_ptr<const FiniteElement> extract_sub_element(const std::vector<uint>& component) const
    
        Extract sub finite element for component

    .. cpp:function:: boost::shared_ptr<const ufc::finite_element> ufc_element() const
    
        Return ufc::finite_element

