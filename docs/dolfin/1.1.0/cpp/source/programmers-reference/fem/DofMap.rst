
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


    .. cpp:function:: DofMap(boost::shared_ptr<const ufc::dofmap> ufc_dofmap, const Mesh& mesh)
    
        Create dof map on mesh (data is not shared)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: DofMap(boost::shared_ptr<const ufc::dofmap> ufc_dofmap, boost::shared_ptr<const Restriction> restriction)
    
        Create restricted dof map on mesh (data is not shared)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            restriction (:cpp:class:`Restriction`)
                The restriction.


    .. cpp:function:: DofMap(const DofMap& dofmap)
    
        Copy constructor
        
        *Arguments*
            dofmap (:cpp:class:`DofMap`)
                The object to be copied.


    .. cpp:function:: bool is_view() const
    
        True iff dof map is a view into another map
        
        *Returns*
            bool
                True if the dof map is a sub-dof map (a view into
                another map).


    .. cpp:function:: bool is_restricted() const
    
        True iff dof map is restricted
        
        *Returns*
            bool
                True iff dof map is restricted


    .. cpp:function:: bool needs_mesh_entities(std::size_t d) const
    
        Return true iff mesh entities of topological dimension d are
        needed
        
        *Arguments*
            d (std::size_t)
                Topological dimension.
        
        *Returns*
            bool
                True if the mesh entities are needed.


    .. cpp:function:: std::size_t global_dimension() const
    
        Return the dimension of the global finite element function
        space
        
        *Returns*
            std::size_t
                The dimension of the global finite element function space.


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


    .. cpp:function:: std::size_t geometric_dimension() const
    
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


    .. cpp:function:: const boost::unordered_map<std::size_t, std::size_t>& off_process_owner() const
    
        Return map from nonlocal dofs that appear in local dof map to
        owning process
        
        *Returns*
            boost::unordered_map<std::size_t, std::size_t>
                The map from non-local dofs.


    .. cpp:function:: const boost::unordered_map<std::size_t, std::vector<std::size_t> >& shared_dofs() const
    
        Return map from all shared dofs to the processes (not including the current
        process) that share it.
        
        *Returns*
            boost::unordered_map<std::size_t, std::vector<std::size_t> >
                The map from dofs to list of processes


    .. cpp:function:: const std::set<std::size_t>& neighbours() const
    
        Return set of all neighbouring processes.
        
        *Returns*
            std::set<std::size_t>
                The set of processes


    .. cpp:function:: const std::vector<dolfin::la_index>& cell_dofs(std::size_t cell_index) const
    
        Local-to-global mapping of dofs on a cell
        
        *Arguments*
            cell_index (std::size_t)
                The cell index.
        
        *Returns*
            std::vector<std::size_t>
                Local-to-global mapping of dofs.


    .. cpp:function:: void tabulate_facet_dofs(unsigned int* dofs, std::size_t local_facet) const
    
        Tabulate local-local facet dofs
        
        *Arguments*
            dofs (std::size_t)
                Degrees of freedom.
            local_facet (std::size_t)
                The local facet.


    .. cpp:function:: void tabulate_coordinates(boost::multi_array<double, 2>& coordinates, const ufc::cell& ufc_cell) const
    
        Tabulate the coordinates of all dofs on a cell (UFC cell
        version)
        
        *Arguments*
            coordinates (boost::multi_array<double, 2>)
                The coordinates of all dofs on a cell.
            ufc_cell (ufc::cell)
                The cell.


    .. cpp:function:: void tabulate_coordinates(boost::multi_array<double, 2>& coordinates, const Cell& cell) const
    
        Tabulate the coordinates of all dofs on a cell (DOLFIN cell
        version)
        
        *Arguments*
            coordinates (boost::multi_array<double, 2>)
                The coordinates of all dofs on a cell.
            cell (:cpp:class:`Cell`)
                The cell.


    .. cpp:function:: boost::shared_ptr<GenericDofMap> copy() const
    
        Create a copy of the dof map
        
        *Returns*
            DofMap
                The Dofmap copy.


    .. cpp:function:: boost::shared_ptr<GenericDofMap> build(const Mesh& new_mesh) const
    
        Create a copy of the dof map
        
        *Arguments*
            new_mesh (:cpp:class:`Mesh`)
                The new mesh to build the dof map on.
        
        *Returns*
            DofMap
                The new Dofmap copy.


    .. cpp:function:: DofMap* extract_sub_dofmap(const std::vector<std::size_t>& component, const Mesh& mesh) const
    
        Extract subdofmap component
        
        *Arguments*
            component (std::vector<std::size_t>)
                The component.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The subdofmap component.


    .. cpp:function:: DofMap* collapse(boost::unordered_map<std::size_t, std::size_t>& collapsed_map, const Mesh& mesh) const
    
        Create a "collapsed" dofmap (collapses a sub-dofmap)
        
        *Arguments*
            collapsed_map (boost::unordered_map<std::size_t, std::size_t>)
                The "collapsed" map.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The collapsed dofmap.


    .. cpp:function:: void set(GenericVector& x, double value) const
    
        Set dof entries in vector to a specified value. Parallel layout
        of vector must be consistent with dof map range.
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector to set.
            value (double)
                The value to set.


    .. cpp:function:: void set_x(GenericVector& x, double value, std::size_t component, const Mesh& mesh) const
    
        Set dof entries in vector to the x[i] coordinate of the dof
        spatial coordinate. Parallel layout of vector must be consistent
        with dof map range.
        
        *Arguments*
            vector (:cpp:class:`GenericVector`)
                The vector to set.
            value (double)
                The value to multiply to coordinate by.
            component (std::size_t)
                The coordinate index.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: boost::unordered_set<std::size_t> dofs() const
    
        Return the set of dof indices
        
        *Returns*
            boost::unordered_set<dolfin::std::size_t>
                The set of dof indices.


    .. cpp:function:: const std::vector<std::vector<dolfin::la_index> >& data() const
    
        Return the underlying dof map data. Intended for internal library
        use only.
        
        *Returns*
            std::vector<std::vector<dolfin::std::size_t> >
                The local-to-global map for each cell.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


