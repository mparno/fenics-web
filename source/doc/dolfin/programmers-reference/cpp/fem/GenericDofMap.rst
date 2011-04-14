.. Documentation for the header file dolfin/fem/GenericDofMap.h

.. _programmers_reference_cpp_fem_genericdofmap:

GenericDofMap.h
===============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: GenericDofMap

    *Parent class*
    
        * :cpp:class:`Variable`
        
    This class provides a generic interface for dof maps

    .. cpp:function:: bool is_view() const = 0
    
        True if dof map is a view into another map (is a sub-dofmap)

    .. cpp:function:: bool needs_mesh_entities(unsigned int d) const = 0
    
        Return true iff mesh entities of topological dimension d are needed

    .. cpp:function:: unsigned int global_dimension() const = 0
    
        Return the dimension of the global finite element function space

    .. cpp:function:: unsigned int cell_dimension(uint index) const = 0
    
        Return the dimension of the local finite element function space on a
        cell

    .. cpp:function:: unsigned int max_cell_dimension() const = 0
    
        Return the maximum dimension of the local finite element function space

    .. cpp:function:: unsigned int num_facet_dofs() const = 0
    
        Return number of facet dofs

    .. cpp:function:: std::pair<unsigned int, unsigned int> ownership_range() const = 0
    
        Return the ownership range (dofs in this range are owned by this process)

    .. cpp:function:: const boost::unordered_map<unsigned int, unsigned int>& off_process_owner() const = 0
    
        Return map from nonlocal-dofs (that appear in local dof map) to owning process

    .. cpp:function:: const std::vector<unsigned int>& cell_dofs(uint cell_index) const = 0
    
        Local-to-global mapping of dofs on a cell

    .. cpp:function:: void tabulate_dofs(uint* dofs, const Cell& cell) const = 0
    
        Tabulate the local-to-global mapping of dofs on a cell

    .. cpp:function:: void tabulate_facet_dofs(uint* dofs, uint local_facet) const = 0
    
        Tabulate local-local facet dofs

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const ufc::cell& ufc_cell) const = 0
    
        Tabulate the coordinates of all dofs on a cell (UFC cell version)

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const Cell& cell) const = 0
    
        Tabulate the coordinates of all dofs on a cell (DOLFIN cell version)

    .. cpp:function:: GenericDofMap* copy(const Mesh& mesh) const = 0
    
        Create a copy of the dof map

    .. cpp:function:: GenericDofMap* extract_sub_dofmap(const std::vector<uint>& component, const Mesh& mesh) const = 0
    
        Extract sub dofmap component

    .. cpp:function:: GenericDofMap* collapse(boost::unordered_map<uint, uint>& collapsed_map, const Mesh& mesh) const = 0
    
        Create a "collapsed" a dofmap (collapses from a sub-dofmap view)

    .. cpp:function:: boost::unordered_set<dolfin::uint> dofs() const = 0
    
        Return the set of dof indices

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

