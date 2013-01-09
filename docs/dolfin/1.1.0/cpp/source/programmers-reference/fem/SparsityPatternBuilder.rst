
.. Documentation for the header file dolfin/fem/SparsityPatternBuilder.h

.. _programmers_reference_cpp_fem_sparsitypatternbuilder:

SparsityPatternBuilder.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SparsityPatternBuilder

    This class provides functions to compute the sparsity pattern.


    .. cpp:function:: static void build(GenericSparsityPattern& sparsity_pattern, const Mesh& mesh, const std::vector<const GenericDofMap*> dofmaps, const std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t> > >& master_slave_dofs, bool cells, bool interior_facets, bool exterior_facets, bool diagonal)
    
        Build sparsity pattern for assembly of given form


