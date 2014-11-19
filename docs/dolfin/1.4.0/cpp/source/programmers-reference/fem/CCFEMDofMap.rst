
.. Documentation for the header file dolfin/fem/CCFEMDofMap.h

.. _programmers_reference_cpp_fem_ccfemdofmap:

CCFEMDofMap.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMDofMap

    *Parent class(es)*
    
        * :cpp:class:`GenericDofMap`
        
    This class handles the mapping of degrees of freedom for CCFEM
    function spaces.


    .. cpp:function:: CCFEMDofMap()
    
        Constructor


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number dofmaps (parts) of the CCFEM dofmap
        
        *Returns*
            std::size_t
                The number of dofmaps (parts) of the CCFEM dofmap


    .. cpp:function:: std::shared_ptr<const GenericDofMap> part(std::size_t i) const
    
        Return dofmap (part) number i
        
        *Returns*
            :cpp:class:`GenericDofMap`
                Dofmap (part) number i


    .. cpp:function:: void set_current_part(std::size_t part) const
    
        Set current part. This will make the CCFEM dofmap act as a
        dofmap for the part of the CCFEM function space defined on the
        current part (mesh).
        
        *Arguments*
            part (std::size_t)
                The number of the part.


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


    .. cpp:function:: void build(const CCFEMFunctionSpace& function_space)
    
        Build CCFEM dofmap


    .. cpp:function:: void clear()
    
        Clear CCFEM dofmap


    .. cpp:function:: bool is_view() const
    
        True if dof map is a view into another map (is a sub-dofmap)


    .. cpp:function:: std::size_t global_dimension() const
    
        Return the dimension of the global finite element function
        space


    .. cpp:function:: std::size_t cell_dimension(std::size_t index) const
    
        Return the dimension of the local finite element function
        space on a cell


    .. cpp:function:: std::size_t max_cell_dimension() const
    
        Return the maximum dimension of the local finite element
        function space


    .. cpp:function:: std::size_t num_entity_dofs(std::size_t dim) const
    
        Return the number of dofs for a given entity dimension


    .. cpp:function:: std::size_t num_facet_dofs() const
    
        Return number of facet dofs


    .. cpp:function:: std::shared_ptr<const Restriction> restriction() const
    
        Restriction if any. If the dofmap is not restricted, a null
        pointer is returned.


    .. cpp:function:: std::pair<std::size_t, std::size_t> ownership_range() const
    
        Return the ownership range (dofs in this range are owned by
        this process)


    .. cpp:function:: const boost::unordered_map<std::size_t, unsigned int>& off_process_owner() const
    
        Return map from nonlocal-dofs (that appear in local dof map)
        to owning process


    .. cpp:function:: const std::vector<dolfin::la_index>& cell_dofs(std::size_t cell_index) const
    
        Local-to-global mapping of dofs on a cell


    .. cpp:function:: void tabulate_facet_dofs(std::vector<std::size_t>& dofs, std::size_t local_facet) const
    
        Tabulate local-local facet dofs


    .. cpp:function:: void tabulate_entity_dofs(std::vector<std::size_t>& dofs, std::size_t dim, std::size_t local_entity) const
    
        Tabulate the local-to-local mapping of dofs on entity
        (dim, local_entity)


    .. cpp:function:: void tabulate_coordinates(boost::multi_array<double, 2>& coordinates, const std::vector<double>& vertex_coordinates, const Cell& cell) const
    
        Tabulate the coordinates of all dofs on a cell (UFC cell version)


    .. cpp:function:: std::vector<double> tabulate_all_coordinates(const Mesh& mesh) const
    
        Tabulate the coordinates of all dofs owned by this
        process. This function is typically used by preconditioners
        that require the spatial coordinates of dofs, for example
        for re-partitioning or nullspace computations. The format for
        the return vector is [x0, y0, z0, x1, y1, z1, . . .].


    .. cpp:function:: std::shared_ptr<GenericDofMap> copy() const
    
        Create a copy of the dof map


    .. cpp:function:: std::shared_ptr<GenericDofMap> create(const Mesh& new_mesh) const
    
        Create a new dof map on new mesh


    .. cpp:function:: std::shared_ptr<GenericDofMap> extract_sub_dofmap(const std::vector<std::size_t>& component, const Mesh& mesh) const
    
        Extract sub dofmap component


    .. cpp:function:: std::shared_ptr<GenericDofMap> collapse(boost::unordered_map<std::size_t, std::size_t>& collapsed_map, const Mesh& mesh) const
    
        Create a "collapsed" a dofmap (collapses from a sub-dofmap view)


    .. cpp:function:: std::vector<dolfin::la_index> dofs() const
    
        Return list of global dof indices on this process


    .. cpp:function:: void set(GenericVector& x, double value) const
    
        Set dof entries in vector to a specified value. Parallel
        layout of vector must be consistent with dof map range. This
        function is typically used to construct the null space of a
        matrix operator


    .. cpp:function:: void set_x(GenericVector& x, double value, std::size_t component, const Mesh& mesh) const
    
        Set dof entries in vector to the value*x[i], where x[i] is the
        spatial coordinate of the dof. Parallel layout of vector must
        be consistent with dof map range. This function is typically
        used to construct the null space of a matrix operator, e.g. rigid
        body rotations.


    .. cpp:function:: const boost::unordered_map<std::size_t, std::vector<unsigned int> >& shared_dofs() const
    
        Return map from shared dofs to the processes (not including
        the current process) that share it.


    .. cpp:function:: const std::set<std::size_t>& neighbours() const
    
        Return set of processes that share dofs with the this process


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


