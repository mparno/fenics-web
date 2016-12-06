
.. Documentation for the header file dolfin/fem/DofMap.h

.. _programmers_reference_cpp_fem_dofmap:

DofMap.h
========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DofMap

    *Parent class(es)*
    
        * :cpp:class:`GenericDofMap`
        
    This class handles the mapping of degrees of freedom. It builds
    a dof map based on a ufc::dofmap on a specific mesh. It will
    reorder the dofs when running in parallel. Sub-dofmaps, both
    views and copies, are supported.


    .. cpp:function:: DofMap(std::shared_ptr<const ufc::dofmap> ufc_dofmap, const Mesh& mesh)
    
        Create dof map on mesh (mesh is not stored)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: DofMap(std::shared_ptr<const ufc::dofmap> ufc_dofmap, const Mesh& mesh, std::shared_ptr<const SubDomain> constrained_domain)
    
        Create a periodic dof map on mesh (mesh is not stored)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            mesh (:cpp:class:`Mesh`)
                The mesh.
            constrained_boundary (:cpp:class:`SubDomain`)
                The subdomain marking the constrained (tied) boundaries.


    .. cpp:function:: bool is_view() const
    
        True iff dof map is a view into another map
        
        *Returns*
            bool
                True if the dof map is a sub-dof map (a view into
                another map).


    .. cpp:function:: std::size_t global_dimension() const
    
        Return the dimension of the global finite element function
        space
        
        *Returns*
            std::size_t
                The dimension of the global finite element function space.


    .. cpp:function:: std::size_t num_element_dofs(std::size_t cell_index) const
    
        Return the dimension of the local finite element function
        space on a cell
        
        *Arguments*
            cell_index (std::size_t)
                Index of cell
        
        *Returns*
            std::size_t
                Dimension of the local finite element function space.


    .. cpp:function:: std::size_t max_element_dofs() const
    
        Return the maximum dimension of the local finite element
        function space
        
        *Returns*
            std::size_t
                Maximum dimension of the local finite element function
                space.


    .. cpp:function:: std::size_t num_entity_dofs(std::size_t dim) const
    
        Return the number of dofs for a given entity dimension
        
        *Arguments*
            dim (std::size_t)
                Entity dimension
        
        *Returns*
            std::size_t
                Number of dofs associated with given entity dimension


    .. cpp:function:: std::size_t num_facet_dofs() const
    
        Return number of facet dofs
        
        *Returns*
            std::size_t
                The number of facet dofs.


    .. cpp:function:: std::pair<std::size_t, std::size_t> ownership_range() const
    
        Return the ownership range (dofs in this range are owned by
        this process)
        
        *Returns*
            std::pair<std::size_t, std::size_t>
                The ownership range.


    .. cpp:function:: const std::vector<int>& off_process_owner() const
    
        Return map from nonlocal dofs that appear in local dof map to
        owning process
        
        *Returns*
            std::vector<unsigned int>
                The map from non-local dofs.


    .. cpp:function:: const std::unordered_map<int, std::vector<int>>& shared_nodes() const
    
        Return map from all shared nodes to the sharing processes (not
        including the current process) that share it.
        
        *Returns*
            std::unordered_map<std::size_t, std::vector<unsigned int>>
                The map from dofs to list of processes


    .. cpp:function:: const std::set<int>& neighbours() const
    
        Return set of processes that share dofs with this process
        
        *Returns*
            std::set<int>
                The set of processes


    .. cpp:function:: void clear_sub_map_data()
    
        Clear any data required to build sub-dofmaps (this is to
        reduce memory use)


    .. cpp:function:: ArrayView<const dolfin::la_index> cell_dofs(std::size_t cell_index) const
    
        Local-to-global mapping of dofs on a cell
        
        *Arguments*
            cell_index (std::size_t)
                The cell index.
        
        *Returns*
            ArrayView<const dolfin::la_index>
                Local-to-global mapping of dofs.


    .. cpp:function:: void tabulate_facet_dofs(std::vector<std::size_t>& dofs, std::size_t local_facet) const
    
        Tabulate local-local facet dofs
        
        *Arguments*
            dofs (std::size_t)
                Degrees of freedom.
            local_facet (std::size_t)
                The local facet.


    .. cpp:function:: void tabulate_entity_dofs(std::vector<std::size_t>& dofs, std::size_t dim, std::size_t local_entity) const
    
        Tabulate local-local mapping of dofs on entity (dim, local_entity)
        
        *Arguments*
            dofs (std::size_t)
                Degrees of freedom.
            dim (std::size_t)
                The entity dimension
            local_entity (std::size_t)
                The local entity index


    .. cpp:function:: void tabulate_global_dofs(std::vector<std::size_t>& dofs) const
    
        Tabulate globally supported dofs
        
        *Arguments*
            dofs (std::size_t)
                Degrees of freedom.


    .. cpp:function:: std::shared_ptr<GenericDofMap> copy() const
    
        Create a copy of the dof map
        
        *Returns*
            DofMap
                The Dofmap copy.


    .. cpp:function:: std::shared_ptr<GenericDofMap> create(const Mesh& new_mesh) const
    
        Create a copy of the dof map on a new mesh
        
        *Arguments*
            new_mesh (:cpp:class:`Mesh`)
                The new mesh to create the dof map on.
        
        *Returns*
            DofMap
                The new Dofmap copy.


    .. cpp:function:: std::shared_ptr<GenericDofMap> extract_sub_dofmap(const std::vector<std::size_t>& component, const Mesh& mesh) const
    
        Extract subdofmap component
        
        *Arguments*
            component (std::vector<std::size_t>)
                The component.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The subdofmap component.


    .. cpp:function:: std::shared_ptr<GenericDofMap> collapse(std::unordered_map<std::size_t, std::size_t>& collapsed_map, const Mesh& mesh) const
    
        Create a "collapsed" dofmap (collapses a sub-dofmap)
        
        *Arguments*
            collapsed_map (std::unordered_map<std::size_t, std::size_t>)
                The "collapsed" map.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The collapsed dofmap.


    .. cpp:function:: std::vector<dolfin::la_index> dofs(const Mesh& mesh, std::size_t dim) const
    
        Return list of dof indices on this process that belong to mesh
        entities of dimension dim


    .. cpp:function:: void set(GenericVector& x, double value) const
    
        Set dof entries in vector to a specified value. Parallel layout
        of vector must be consistent with dof map range. This
        function is typically used to construct the null space of a
        matrix operator.
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector to set.
            value (double)
                The value to set.


    .. cpp:function:: std::shared_ptr<const IndexMap> index_map() const
    
        Return the map (const access)


    .. cpp:function:: int block_size() const
    
        Return the block size for dof maps with components, typically
        used for vector valued functions.


    .. cpp:function:: void tabulate_local_to_global_dofs(std::vector<std::size_t>& local_to_global_map) const
    
        Compute the map from local (this process) dof indices to
        global dof indices.
        
        *Arguments*
            local_to_global_map (_std::vector<std::size_t>_)
                The local-to-global map to fill.


    .. cpp:function:: std::size_t local_to_global_index(int local_index) const
    
        Return global dof index for a given local (process) dof index
        
        *Arguments*
            local_index (int)
                The local local index.
        
        *Returns*
            std::size_t
                The global dof index.


    .. cpp:function:: const std::vector<std::size_t>& local_to_global_unowned() const
    
        Return indices of dofs which are owned by other processes


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


