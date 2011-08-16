
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


    .. cpp:function:: DofMap(boost::shared_ptr<const ufc::dofmap> ufc_dofmap, Mesh& mesh)
    
        Create dof map on mesh (data is not shared)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: DofMap(boost::shared_ptr<const ufc::dofmap> ufc_dofmap, const Mesh& mesh)
    
        Create dof map on mesh ((data is not shared), const mesh
        version)
        
        *Arguments*
            ufc_dofmap (ufc::dofmap)
                The ufc::dofmap.
            mesh (:cpp:class:`Mesh`)
                The mesh.


    .. cpp:function:: DofMap(const DofMap& dofmap)
    
        Copy constructor
        
        *Arguments*
            dofmap (:cpp:class:`DofMap`)
                The object to be copied.


    .. cpp:function:: DofMap(const DofMap& parent_dofmap, const std::vector<uint>& component, const Mesh& mesh, bool distributed)
    
        Create a sub-dofmap (a view) from parent_dofmap


    .. cpp:function:: DofMap(boost::unordered_map<uint, uint>& collapsed_map, const DofMap& dofmap_view, const Mesh& mesh, bool distributed)
    
        Create a collapsed dofmap from parent_dofmap


    .. cpp:function:: bool is_view() const
    
        True if dof map is a view into another map
        
        *Returns*
            bool
                True if the dof map is a sub-dof map (a view into
                another map).


    .. cpp:function:: bool needs_mesh_entities(unsigned int d) const
    
        Return true iff mesh entities of topological dimension d are
        needed
        
        *Arguments*
            d (unsigned int)
                Topological dimension.
        
        *Returns*
            bool
                True if the mesh entities are needed.


    .. cpp:function:: unsigned int global_dimension() const
    
        Return the dimension of the global finite element function
        space
        
        *Returns*
            unsigned int
                The dimension of the global finite element function space.


    .. cpp:function:: unsigned int cell_dimension(uint cell_index) const
    
        Return the dimension of the local finite element function
        space on a cell
        
        *Arguments*
            cell_index (uint)
                Index of cell
        
        *Returns*
            unsigned int
                Dimension of the local finite element function space.


    .. cpp:function:: unsigned int max_cell_dimension() const
    
        Return the maximum dimension of the local finite element
        function space
        
        *Returns*
            unsigned int
                Maximum dimension of the local finite element function
                space.


    .. cpp:function:: unsigned int geometric_dimension() const
    
        Return the geometric dimension of the coordinates this dof map
        provides
        
        *Returns*
            unsigned int
                The geometric dimension.


    .. cpp:function:: unsigned int num_facet_dofs() const
    
        Return number of facet dofs
        
        *Returns*
            unsigned int
                The number of facet dofs.


    .. cpp:function:: std::pair<unsigned int, unsigned int> ownership_range() const
    
        Return the ownership range (dofs in this range are owned by
        this process)
        
        *Returns*
            std::pair<unsigned int, unsigned int>
                The ownership range.


    .. cpp:function:: const boost::unordered_map<unsigned int, unsigned int>& off_process_owner() const
    
        Return map from nonlocal dofs that appear in local dof map to
        owning process
        
        *Returns*
            boost::unordered_map<unsigned int, unsigned int>
                The map from non-local dofs.


    .. cpp:function:: const std::vector<uint>& cell_dofs(uint cell_index) const
    
        Local-to-global mapping of dofs on a cell
        
        *Arguments*
            cell_index (uint)
                The cell index.
        
        *Returns*
            std::vector<uint>
                Local-to-global mapping of dofs.


    .. cpp:function:: void tabulate_dofs(uint* dofs, const Cell& cell) const
    
        Tabulate the local-to-global mapping of dofs on a cell
        
        *Arguments*
            dofs (uint)
                Degrees of freedom on a cell.
            cell (:cpp:class:`Cell`)
                The cell.


    .. cpp:function:: void tabulate_facet_dofs(uint* dofs, uint local_facet) const
    
        Tabulate local-local facet dofs
        
        *Arguments*
            dofs (uint)
                Degrees of freedom.
            local_facet (uint)
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


    .. cpp:function:: DofMap* copy(const Mesh& mesh) const
    
        Create a copy of the dof map
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The object to be copied.


    .. cpp:function:: DofMap* extract_sub_dofmap(const std::vector<uint>& component, const Mesh& mesh) const
    
        Extract subdofmap component
        
        *Arguments*
            component (std::vector<uint>)
                The component.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The subdofmap component.


    .. cpp:function:: DofMap* collapse(boost::unordered_map<uint, uint>& collapsed_map, const Mesh& mesh) const
    
        Create a "collapsed" dofmap (collapses a sub-dofmap)
        
        *Arguments*
            collapsed_map (boost::unordered_map<uint, uint>)
                The "collapsed" map.
            mesh (:cpp:class:`Mesh`)
                The mesh.
        
        *Returns*
            DofMap
                The collapsed dofmap.


    .. cpp:function:: boost::unordered_set<dolfin::uint> dofs() const
    
        Return the set of dof indices
        
        *Returns*
            boost::unordered_set<dolfin::uint>
                The set of dof indices.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


