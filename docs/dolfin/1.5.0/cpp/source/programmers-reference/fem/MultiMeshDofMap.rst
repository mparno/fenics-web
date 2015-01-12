
.. Documentation for the header file dolfin/fem/MultiMeshDofMap.h

.. _programmers_reference_cpp_fem_multimeshdofmap:

MultiMeshDofMap.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshDofMap

    This class handles the mapping of degrees of freedom for MultiMesh
    function spaces.


    .. cpp:function:: MultiMeshDofMap()
    
        Constructor


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number dofmaps (parts) of the MultiMesh dofmap
        
        *Returns*
            std::size_t
                The number of dofmaps (parts) of the MultiMesh dofmap


    .. cpp:function:: std::shared_ptr<const GenericDofMap> part(std::size_t i) const
    
        Return dofmap (part) number i
        
        *Returns*
            :cpp:class:`GenericDofMap`
                Dofmap (part) number i


    .. cpp:function:: void add(std::shared_ptr<const GenericDofMap> dofmap)
    
        Add dofmap (shared pointer version)
        
        *Arguments*
            dofmap (:cpp:class:`GenericDofMap`)
                The dofmap.


    .. cpp:function:: void add(const GenericDofMap& dofmap)
    
        Add dofmap (reference version)
        
        *Arguments*
            dofmap (:cpp:class:`DofMap`)
                The dofmap.


    .. cpp:function:: void build(const MultiMeshFunctionSpace& function_space, const std::vector<dolfin::la_index>& offsets)
    
        Build MultiMesh dofmap


    .. cpp:function:: void clear()
    
        Clear MultiMesh dofmap


    .. cpp:function:: std::size_t global_dimension() const
    
        Return the dimension of the global finite element function
        space


    .. cpp:function:: std::pair<std::size_t, std::size_t> ownership_range() const
    
        Return the ownership range (dofs in this range are owned by
        this process)


    .. cpp:function:: const std::vector<int>& off_process_owner() const
    
        Return map from nonlocal-dofs (that appear in local dof map)
        to owning process


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


