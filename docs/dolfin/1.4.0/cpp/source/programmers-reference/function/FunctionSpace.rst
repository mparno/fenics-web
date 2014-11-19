
.. Documentation for the header file dolfin/function/FunctionSpace.h

.. _programmers_reference_cpp_function_functionspace:

FunctionSpace.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FunctionSpace

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
        * :cpp:class:`Hierarchical<FunctionSpace>`
        
    This class represents a finite element function space defined by
    a mesh, a finite element, and a local-to-global mapping of the
    degrees of freedom (dofmap).


    .. cpp:function:: FunctionSpace(std::shared_ptr<const Mesh> mesh, std::shared_ptr<const FiniteElement> element, std::shared_ptr<const GenericDofMap> dofmap)
    
        Create function space for given mesh, element and dofmap
        (shared data)
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            element (:cpp:class:`FiniteElement`)
                The element.
            dofmap (:cpp:class:`GenericDofMap`)
                The dofmap.


    .. cpp:function:: explicit FunctionSpace(std::shared_ptr<const Mesh> mesh)
    
        Create empty function space for later initialization. This
        constructor is intended for use by any sub-classes which need
        to construct objects before the initialisation of the base
        class. Data can be attached to the base class using
        FunctionSpace::attach(...).
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: FunctionSpace(const FunctionSpace& V)
    
        Copy constructor
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The object to be copied.


    .. cpp:function:: void attach(std::shared_ptr<const FiniteElement> element, std::shared_ptr<const GenericDofMap> dofmap)
    
        Attach data to an empty function space
        
        *Arguments*
            element (:cpp:class:`FiniteElement`)
                The element.
            dofmap (:cpp:class:`GenericDofMap`)
                The dofmap.


    .. cpp:function:: const FunctionSpace& operator= (const FunctionSpace& V)
    
        Assignment operator
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                Another function space.


    .. cpp:function:: bool operator== (const FunctionSpace& V) const
    
        Equality operator
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                Another function space.


    .. cpp:function:: bool operator!= (const FunctionSpace& V) const
    
        Unequality operator
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                Another function space.


    .. cpp:function:: std::shared_ptr<const Mesh> mesh() const
    
        Return mesh
        
        *Returns*
            :cpp:class:`Mesh`
                The mesh.


    .. cpp:function:: std::shared_ptr<const FiniteElement> element() const
    
        Return finite element
        
        *Returns*
            :cpp:class:`FiniteElement`
                The finite element.


    .. cpp:function:: std::shared_ptr<const GenericDofMap> dofmap() const
    
        Return dofmap
        
        *Returns*
            :cpp:class:`GenericDofMap`
                The dofmap.


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of function space
        
        *Returns*
            std::size_t
                The dimension of the function space.


    .. cpp:function:: void interpolate(GenericVector& expansion_coefficients, const GenericFunction& v) const
    
        Interpolate function v into function space, returning the
        vector of expansion coefficients
        
        *Arguments*
            expansion_coefficients (:cpp:class:`GenericVector`)
                The expansion coefficients.
            v (:cpp:class:`GenericFunction`)
                The function to be interpolated.


    .. cpp:function:: std::shared_ptr<FunctionSpace> operator[] (std::size_t i) const
    
        Extract subspace for component
        
        *Arguments*
            i (std::size_t)
                Index of the subspace.
        *Returns*
            :cpp:class:`FunctionSpace`
                The subspace.


    .. cpp:function:: std::shared_ptr<FunctionSpace> extract_sub_space(const std::vector<std::size_t>& component) const
    
        Extract subspace for component
        
        *Arguments*
            component (std::vector<std::size_t>)
                The component.
        
        *Returns*
            :cpp:class:`FunctionSpace`
                The subspace.


    .. cpp:function:: std::shared_ptr<FunctionSpace> collapse() const
    
        Collapse a subspace and return a new function space
        
        *Returns*
            :cpp:class:`FunctionSpace`
                The new function space.


    .. cpp:function:: std::shared_ptr<FunctionSpace> collapse(boost::unordered_map<std::size_t, std::size_t>& collapsed_dofs) const
    
        Collapse a subspace and return a new function space and a map
        from new to old dofs
        
        *Arguments*
            collapsed_dofs (boost::unordered_map<std::size_t, std::size_t>)
                The map from new to old dofs.
        
        *Returns*
            :cpp:class:`FunctionSpace`
              The new function space.


    .. cpp:function:: bool has_cell(const Cell& cell) const
    
        Check if function space has given cell
        
        *Arguments*
            cell (:cpp:class:`Cell`)
                The cell.
        
        *Returns*
            bool
                True if the function space has the given cell.


    .. cpp:function:: bool has_element(const FiniteElement& element) const
    
        Check if function space has given element
        
        *Arguments*
            element (:cpp:class:`FiniteElement`)
                The finite element.
        
        *Returns*
            bool
                True if the function space has the given element.


    .. cpp:function:: std::vector<std::size_t> component() const
    
        Return component
        
        *Returns*
            std::vector<std::size_t>
                The component (relative to superspace).


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


    .. cpp:function:: void print_dofmap() const
    
        Print dofmap (useful for debugging)


