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

    .. cpp:function:: GenericDofMap* collapse(std::map<uint, uint>& collapsed_map, const Mesh& dolfin_mesh) const = 0
    
        "Collapse" a sub dofmap

    .. cpp:function:: GenericDofMap* extract_sub_dofmap(const std::vector<uint>& component, const Mesh& dolfin_mesh) const = 0
    
        Extract sub dofmap component

    .. cpp:function:: Set<dolfin::uint> dofs(const Mesh& mesh, bool sort = false) const = 0
    
        Return the set of dof indices

    .. cpp:function:: bool needs_mesh_entities(unsigned int d) const = 0
    
        Return true iff mesh entities of topological dimension d are needed

    .. cpp:function:: std::string signature() const = 0
    
        Return a string identifying the dof map

    .. cpp:function:: std::string str(bool verbose) const = 0
    
        Return informal string representation (pretty-print)

    .. cpp:function:: unsigned int global_dimension() const = 0
    
        Return the dimension of the global finite element function space

    .. cpp:function:: unsigned int local_dimension(const ufc::cell& cell) const = 0
    
        Return the dimension of the local finite element function space on a
        cell

    .. cpp:function:: unsigned int max_local_dimension() const = 0
    
        Return the maximum dimension of the local finite element function space

    .. cpp:function:: unsigned int num_facet_dofs() const = 0
    
        Return number of facet dofs

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const Cell& cell) const = 0
    
        Tabulate the coordinates of all dofs on a cell (DOLFIN cell version)

    .. cpp:function:: void tabulate_coordinates(double** coordinates, const ufc::cell& ufc_cell) const = 0
    
        Tabulate the coordinates of all dofs on a cell (UFC cell version)

    .. cpp:function:: void tabulate_dofs(uint* dofs, const Cell& cell) const = 0
    
        Tabulate the local-to-global mapping of dofs on a cell
        (DOLFIN cell version)

    .. cpp:function:: void tabulate_dofs(uint* dofs, const ufc::cell& ufc_cell, uint cell_index) const = 0
    
        Tabulate the local-to-global mapping of dofs on a cell
        (UFC cell version)

    .. cpp:function:: void tabulate_facet_dofs(uint* dofs, uint local_facet) const = 0
    
        Tabulate local-local facet dofs

