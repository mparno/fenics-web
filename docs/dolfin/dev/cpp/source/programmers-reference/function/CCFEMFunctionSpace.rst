
.. Documentation for the header file dolfin/function/CCFEMFunctionSpace.h

.. _programmers_reference_cpp_function_ccfemfunctionspace:

CCFEMFunctionSpace.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMFunctionSpace

    This class represents a cut and composite finite element
    function space (CCFEM) defined on one or more possibly
    intersecting meshes.
    
    FIXME: Document usage of class with add() followed by build()


    .. cpp:function:: CCFEMFunctionSpace()
    
        Create empty CCFEM function space


    .. cpp:function:: std::size_t dim() const
    
        Return dimension of the CCFEM function space
        
        *Returns*
            std::size_t
                The dimension of the CCFEM function space.


    .. cpp:function:: boost::shared_ptr<const CCFEMDofMap> dofmap() const
    
        Return CCFEM dofmap
        
        *Returns*
            :cpp:class:`CCFEMDofMap`
                The dofmap.


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number function spaces (parts) of the CCFEM function space
        
        *Returns*
            std::size_t
                The number of function spaces (parts) of the CCFEM function space.


    .. cpp:function:: boost::shared_ptr<const FunctionSpace> part(std::size_t i) const
    
        Return function space (part) number i
        
        *Returns*
            :cpp:class:`FunctionSpace`
                Function space (part) number i


    .. cpp:function:: void add(boost::shared_ptr<const FunctionSpace> function_space)
    
        Add function space (shared pointer version)
        
        *Arguments*
            function_space (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: void add(const FunctionSpace& function_space)
    
        Add function space (reference version)
        
        *Arguments*
            function_space (:cpp:class:`FunctionSpace`)
                The function space.


    .. cpp:function:: void build()
    
        Build CCFEM function space


    .. cpp:function:: void clear()
    
        Clear CCFEM function space


