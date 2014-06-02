
.. Documentation for the header file dolfin/fem/SparsityPatternBuilder.h

.. _programmers_reference_cpp_fem_sparsitypatternbuilder:

SparsityPatternBuilder.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SparsityPatternBuilder

    This class provides functions to compute the sparsity pattern.


    .. cpp:function:: static void build(GenericSparsityPattern& sparsity_pattern, const Mesh& mesh, const std::vector<const GenericDofMap*> dofmaps, bool cells, bool interior_facets, bool exterior_facets, bool diagonal, bool init=true, bool finalize=true)
    
        Build sparsity pattern for assembly of given form


    .. cpp:function:: static void build_ccfem(GenericSparsityPattern& sparsity_pattern, const CCFEMForm& form)
    
        Build sparsity pattern for assembly of given CCFEM form


