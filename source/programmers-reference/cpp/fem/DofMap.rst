.. Documentation for the header file dolfin/fem/DofMap.h

.. _programmers_reference_cpp_fem_dofmap:

DofMap.h
========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: DofMap

    *Parent class*
    
        * :cpp:class:`GenericDofMap`
        
        This class handles the mapping of degrees of freedom. It builds
        a dof map based on a ufc::dof_map on a specific mesh. It will
        reorder the dofs when running in parallel.
        
        If ufc_offset != 0, then the dof map provides a view into a
        larger dof map. A dof map which is a view, can be 'collapsed'
        such that the dof indices are contiguous.

    .. cpp:function:: DofMap(boost::shared_ptr<ufc::dof_map> ufc_dofmap,
                             const Mesh& dolfin_mesh)
    
        Create dof map on mesh (const mesh version)

    .. cpp:function:: DofMap(boost::shared_ptr<ufc::dof_map> ufc_dofmap,
                       Mesh& dolfin_mesh)
    
        Create dof map on mesh

    .. cpp:function:: DofMap(boost::shared_ptr<ufc::dof_map> ufc_dofmap, const UFCMesh& ufc_mesh)
    
        Create dof map on mesh with a std::vector dof map

    .. cpp:function:: DofMap* collapse(std::map<uint, uint>& collapsed_map, const Mesh& dolfin_mesh) const
    
        "Collapse" a sub dofmap

    .. cpp:function:: DofMap* extract_sub_dofmap(const std::vector<uint>& component, const Mesh& dolfin_mesh) const
    
        Extract sub dofmap component

    .. cpp:function:: Set<dolfin::uint> dofs(const Mesh& mesh, bool sort = false) const
    
        Return the set of dof indices

    .. cpp:function:: bool needs_mesh_entities(unsigned int d) const
    
        Return true iff mesh entities of topological dimension d are needed

    .. cpp:function:: friend class DofMapBuilder
    
        Friends

    .. cpp:function:: static void init_ufc_dofmap(ufc::dof_map& dofmap,
                                                  const ufc::mesh ufc_mesh,
                                                  const Mesh& dolfin_mesh)
    
        Initialize the UFC dofmap

    .. cpp:function:: std::string signature() const
    
        Return a string identifying the dof map

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: unsigned int global_dimension() const
    
        Return the dimension of the global finite element function space

    .. cpp:function:: unsigned int local_dimension(const ufc::cell& cell) const
    
        Return the dimension of the local finite element function space on a cell

    .. cpp:function:: unsigned int max_local_dimension() const
    
        Return the maximum dimension of the local finite element function space

    .. cpp:function:: unsigned int num_facet_dofs() const
    
        Return number of facet dofs

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const Cell& cell) const
    
        Tabulate the coordinates of all dofs on a cell (DOLFIN cell version)

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const ufc::cell& ufc_cell) const
    
        Tabulate the coordinates of all dofs on a cell (UFC cell version)

    .. cpp:function:: void tabulate_dofs(uint* dofs, const Cell& cell) const
    
        Tabulate the local-to-global mapping of dofs on a cell (DOLFIN cell version)

    .. cpp:function:: void tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell, uint cell_index) const
    
        Tabulate the local-to-global mapping of dofs on a cell (UFC cell version)

    .. cpp:function:: void tabulate_facet_dofs(uint* dofs, uint local_facet) const
    
        Tabulate local-local facet dofs

    .. cpp:function:: ~DofMap()
    
        Destructor

