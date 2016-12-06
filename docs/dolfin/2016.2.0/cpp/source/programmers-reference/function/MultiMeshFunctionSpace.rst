
.. Documentation for the header file dolfin/function/MultiMeshFunctionSpace.h

.. _programmers_reference_cpp_function_multimeshfunctionspace:

MultiMeshFunctionSpace.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshFunctionSpace

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class represents a function space on a multimesh. It may
    may be created from a set of standard function spaces by
    repeatedly calling add(), followed by a call to build(). Note
    that a multimesh function space is not useful and its data
    structures are empty until build() has been called.


    .. cpp:function:: MultiMeshFunctionSpace(std::shared_ptr<const MultiMesh> multimesh)
    
        Create multimesh function space on multimesh (shared pointer version)


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of the multimesh function space
        
        *Returns*
            std::size_t
                The dimension of the multimesh function space.


    .. cpp:function:: std::shared_ptr<const MultiMesh> multimesh() const
    
        Return multimesh
        
        *Returns*
            :cpp:class:`MultiMesh`
                The multimesh.


    .. cpp:function:: std::shared_ptr<const MultiMeshDofMap> dofmap() const
    
        Return multimesh dofmap
        
        *Returns*
            :cpp:class:`MultiMeshDofMap`
                The dofmap.


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number of function spaces (parts) of the multimesh function space
        
        *Returns*
            std::size_t
                The number of function spaces (parts) of the multimesh function space.


    .. cpp:function:: std::shared_ptr<const FunctionSpace> part(std::size_t i) const
    
        Return function space (part) number i
        
        *Arguments*
            i (std::size_t)
                The part number
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Function space (part) number i


    .. cpp:function:: std::shared_ptr<const FunctionSpace> view(std::size_t i) const
    
        Return view of multimesh function space for part number i.
        This function differs from the part() function in that it does
        not return the original function space for a part, but rather
        a view of the common multimesh function space (dofs global to
        the collection of parts).
        
        *Arguments*
            i (std::size_t)
                The part number
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Function space (part) number i


    .. cpp:function:: void add(std::shared_ptr<const FunctionSpace> function_space)
    
        Add function space
        
        *Arguments*
            function_space (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: void build()
    
        Build multimesh function space


    .. cpp:function:: void build(const std::vector<dolfin::la_index>& offsets)
    
        Build multimesh function space. This function uses offsets
        computed from the full function spaces on each part.


