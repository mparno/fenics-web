.. Documentation for the header file dolfin/function/FunctionSpace.h

.. _programmers_reference_cpp_function_functionspace:

FunctionSpace.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: FunctionSpace

    *Parent class*
    
        * :cpp:class:`Variable`
        
    This class represents a finite element function space defined by
    a mesh, a finite element, and a local-to-global mapping of the
    degrees of freedom (dofmap).

    .. cpp:function:: FunctionSpace(boost::shared_ptr<Mesh> mesh)
    
        Create an empty function space for later intialisation. This
        constructor is intended for use by any sub-classes which need to
        construct objects before the initialisation of the base class. Data can
        be attached to the base class using FunctionSpace::attach(...)

    .. cpp:function:: FunctionSpace(boost::shared_ptr<Mesh> mesh, boost::shared_ptr<const FiniteElement> element, boost::shared_ptr<const GenericDofMap> dofmap)
    
        Create function space for given mesh, element and dofmap (shared data)

    .. cpp:function:: FunctionSpace(boost::shared_ptr<const Mesh> mesh, boost::shared_ptr<const FiniteElement> element, boost::shared_ptr<const GenericDofMap> dofmap)
    
        Create function space for given mesh, element and dofmap (shared data)

    .. cpp:function:: FunctionSpace(const FunctionSpace& V)
    
        Copy constructor

    .. cpp:function:: bool has_cell(const Cell& cell) const
    
        Check if function space has given cell

    .. cpp:function:: bool has_element(const FiniteElement& element) const
    
        Check if function space has given element

    .. cpp:function:: boost::shared_ptr<FunctionSpace> collapse_sub_space(boost::shared_ptr<GenericDofMap> dofmap) const
    
        Return function space with a new dof map

    .. cpp:function:: boost::shared_ptr<FunctionSpace> extract_sub_space(const std::vector<uint>& component) const
    
        Extract sub space for component

    .. cpp:function:: boost::shared_ptr<FunctionSpace> operator[] (uint i) const
    
        Extract sub space for component

    .. cpp:function:: const Array<uint>& component() const
    
        Return component (relative to super space)

    .. cpp:function:: const FiniteElement& element() const
    
        Return finite element

    .. cpp:function:: const FunctionSpace& operator= (const FunctionSpace& V)
    
        Assignment operator

    .. cpp:function:: const GenericDofMap& dofmap() const
    
        Return dofmap

    .. cpp:function:: const Mesh& mesh() const
    
        Return mesh

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: uint dim() const
    
        Return dimension of function space

    .. cpp:function:: void attach(boost::shared_ptr<const FiniteElement> element, boost::shared_ptr<const GenericDofMap> dofmap)
    
        Attach data to an empty FunctionSpace

    .. cpp:function:: void interpolate(GenericVector& expansion_coefficients, const GenericFunction& v) const
    
        Interpolate function v into function space, returning the vector of
        expansion coefficients

    .. cpp:function:: void print_dofmap() const
    
        Print dofmap (useful for debugging)

