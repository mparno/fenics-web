
.. Documentation for the header file dolfin/fem/GenericDofMap.h

.. _programmers_reference_cpp_fem_genericdofmap:

GenericDofMap.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericDofMap

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class provides a generic interface for dof maps


    .. cpp:function:: GenericDofMap()
    
        Constructor


    .. cpp:function:: bool is_view() const = 0
    
        True if dof map is a view into another map (is a sub-dofmap)


    .. cpp:function:: std::size_t global_dimension() const = 0
    
        Return the dimension of the global finite element function
        space


    .. cpp:function:: std::size_t cell_dimension(std::size_t index) const
    
        Return the dimension of the local finite element function
        space on a cell (deprecated API)


    .. cpp:function:: std::size_t max_cell_dimension() const
    
        Return the maximum dimension of the local finite element
        function space (deprecated API)


    .. cpp:function:: std::size_t num_element_dofs(std::size_t index) const = 0
    
        Return the dimension of the local finite element function
        space on a cell


    .. cpp:function:: std::size_t max_element_dofs() const = 0
    
        Return the maximum dimension of the local finite element
        function space


    .. cpp:function:: std::size_t num_entity_dofs(std::size_t dim) const = 0
    
        Return the number of dofs for a given entity dimension


    .. cpp:function:: std::size_t num_facet_dofs() const = 0
    
        Return number of facet dofs


    .. cpp:function:: std::pair<std::size_t, std::size_t> ownership_range() const = 0
    
        Return the ownership range (dofs in this range are owned by
        this process)


    .. cpp:function:: const std::vector<int>& off_process_owner() const = 0
    
        Return map from nonlocal-dofs (that appear in local dof map)
        to owning process


    .. cpp:function:: ArrayView<const dolfin::la_index> cell_dofs(std::size_t cell_index) const = 0
    
        Local-to-global mapping of dofs on a cell


    .. cpp:function:: void tabulate_facet_dofs(std::vector<std::size_t>& dofs, std::size_t local_facet) const = 0
    
        Tabulate local-local facet dofs


    .. cpp:function:: void tabulate_entity_dofs(std::vector<std::size_t>& dofs, std::size_t dim, std::size_t local_entity) const = 0
    
        Tabulate the local-to-local mapping of dofs on entity
        (dim, local_entity)


    .. cpp:function:: void tabulate_global_dofs(std::vector<std::size_t>& dofs) const = 0
    
        Tabulate globally supported dofs


    .. cpp:function:: std::shared_ptr<GenericDofMap> copy() const = 0
    
        Create a copy of the dof map


    .. cpp:function:: std::shared_ptr<GenericDofMap> create(const Mesh& new_mesh) const = 0
    
        Create a new dof map on new mesh


    .. cpp:function:: std::shared_ptr<GenericDofMap> extract_sub_dofmap(const std::vector<std::size_t>& component, const Mesh& mesh) const = 0
    
        Extract sub dofmap component


    .. cpp:function:: std::shared_ptr<GenericDofMap> collapse(std::unordered_map<std::size_t, std::size_t>& collapsed_map, const Mesh& mesh) const = 0
    
        Create a "collapsed" a dofmap (collapses from a sub-dofmap view)


    .. cpp:function:: std::vector<dolfin::la_index> dofs(const Mesh& mesh, std::size_t dim) const = 0
    
        Return list of dof indices on this process that belong to mesh
        entities of dimension dim


    .. cpp:function:: std::vector<dolfin::la_index> dofs() const = 0
    
        Return list of global dof indices on this process


    .. cpp:function:: void set(GenericVector& x, double value) const = 0
    
        Set dof entries in vector to a specified value. Parallel
        layout of vector must be consistent with dof map range. This
        function is typically used to construct the null space of a
        matrix operator


    .. cpp:function:: const std::vector<std::size_t>& local_to_global_unowned() const = 0
    
        Return the map from unowned local dofmap nodes to global dofmap
        nodes. Dofmap node is dof index modulo block size.


    .. cpp:function:: std::shared_ptr<const IndexMap> index_map() const = 0
    
        Index map (const access)


    .. cpp:function:: void tabulate_local_to_global_dofs(std::vector<std::size_t>& local_to_global_map) const = 0
    
        Tabulate map between local (process) and global dof indices


    .. cpp:function:: std::size_t local_to_global_index(int local_index) const = 0
    
        Return global dof index corresponding to a given local index


    .. cpp:function:: const std::unordered_map<int, std::vector<int>>& shared_nodes() const = 0
    
        Return map from shared nodes to the processes (not including
        the current process) that share it.


    .. cpp:function:: const std::set<int>& neighbours() const = 0
    
        Return set of processes that share dofs with the this process


    .. cpp:function:: void clear_sub_map_data() = 0
    
        Clear any data required to build sub-dofmaps (this is to
        reduce memory use)


    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)


    .. cpp:function:: int block_size() const = 0
    
        Get block size


