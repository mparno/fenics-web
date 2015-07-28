
.. Documentation for the header file dolfin/fem/SparsityPatternBuilder.h

.. _programmers_reference_cpp_fem_sparsitypatternbuilder:

SparsityPatternBuilder.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SparsityPatternBuilder

    This class provides functions to compute the sparsity pattern
    based on DOF maps


    .. cpp:function:: static void build(GenericSparsityPattern& sparsity_pattern, const Mesh& mesh, const std::vector<const GenericDofMap*> dofmaps, bool cells, bool interior_facets, bool exterior_facets, bool vertices, bool diagonal, bool init=true, bool finalize=true)
    
        Build sparsity pattern for assembly of given form


    .. cpp:function:: static void build_multimesh_sparsity_pattern (GenericSparsityPattern& sparsity_pattern, const MultiMeshForm& form)
    
        Build sparsity pattern for assembly of given multimesh form


    .. cpp:function:: static void _build_multimesh_sparsity_pattern_interface (GenericSparsityPattern& sparsity_pattern, const MultiMeshForm& form, std::size_t part)
    
        Build sparsity pattern for interface part of multimesh form


