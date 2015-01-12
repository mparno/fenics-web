
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


    .. cpp:function:: std::size_t local_dimension(std::string type) const
    
        Return number of owned, unowned, or all dofs in the dofmap on
        this process
        
        *Arguments*
            type (std::string)
                Either "owned", "unowned", or "all"
        
        *Returns*
            std::size_t
                Number of local dofs.


    .. cpp:function:: std::size_t cell_dimension(std::size_t cell_index) const
    
        Return the dimension of the local finite element function
        space on a cell
        
        *Arguments*
            cell_index (std::size_t)
                Index of cell
        
        *Returns*
            std::size_t
                Dimension of the local finite element function space.


    .. cpp:function:: std::size_t max_cell_dimension() const
    
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


    .. cpp:function:: std::size_t geometric_dimension() const
    
        DEPRECATED
        Return the geometric dimension of the coordinates this dof map
        provides
        
        *Returns*
            std::size_t
                The geometric dimension.


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


    .. cpp:function:: const std::unordered_map<int, std::vector<int> >& shared_nodes() const
    
        Return map from all shared nodes to the sharing processes (not
        including the current process) that share it.
        
        *Returns*
            std::unordered_map<std::size_t, std::vector<unsigned int> >
                The map from dofs to list of processes


    .. cpp:function:: const std::set<int>& neighbours() const
    
        Return set of processes that share dofs with this process
        
        *Returns*
            std::set<int>
                The set of processes


    .. cpp:function:: void clear_sub_map_data()
    
        Clear any data required to build sub-dofmaps (this is to
        reduce memory use)


    .. cpp:function:: const std::vector<dolfin::la_index>& cell_dofs(std::size_t cell_index) const
    
        Local-to-global mapping of dofs on a cell
        
        *Arguments*
            cell_index (std::size_t)
                The cell index.
        
        *Returns*
            std::vector<dolfin::la_index>
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


    .. cpp:function:: void tabulate_coordinates(boost::multi_array<double, 2>& coordinates, const std::vector<double>& vertex_coordinates, const Cell& cell) const
    
        Tabulate the coordinates of all dofs on a cell (UFC cell
        version)
        
        *Arguments*
            coordinates (boost::multi_array<double, 2>)
                The coordinates of all dofs on a cell.
            vertex_coordinates (std::vector<double>)
                The cell vertex coordinates
            cell (Cell)
                The cell.


    .. cpp:function:: std::vector<double> tabulate_all_coordinates(const Mesh& mesh) const
    
        Tabulate the coordinates of all dofs on this process. This
        function is typically used by preconditioners that require the
        spatial coordinates of dofs, for example for re-partitioning or
        nullspace computations.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            std::vector<double>
                The dof coordinates (x0, y0, x1, y1, . . .)


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


    .. cpp:function:: void set_x(GenericVector& x, double value, std::size_t component, const Mesh& mesh) const
    
        Set dof entries in vector to the x[i] coordinate of the dof
        spatial coordinate. Parallel layout of vector must be consistent
        with dof map range This function is typically used to
        construct the null space of a matrix operator, e.g. rigid
        body rotations.
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector to set.
            value (double)
                The value to multiply to coordinate by.
            component (std::size_t)
                The coordinate index.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: const std::vector<std::size_t>& local_to_global_unowned() const
    
        Return the map from unowned local dofmap nodes to global dofmap
        nodes. Dofmap node is dof index modulo block size.
        
        *Returns*
            _std::vector<std::size_t>_
                The unonwed local-to-global node map.


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


    .. cpp:function:: const std::vector<std::vector<dolfin::la_index> >& data() const
    
        Return the underlying dof map data. Intended for internal library
        use only.
        
        *Returns*
            std::vector<std::vector<dolfin::la_index> >
                The local-to-global map for each cell.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


