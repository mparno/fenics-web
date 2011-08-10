
.. Documentation for the header file dolfin/adaptivity/LocalAssembler.h

.. _programmers_reference_cpp_adaptivity_localassembler:

LocalAssembler.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalAssembler

    


    .. cpp:function:: static void assemble(arma::mat& A, UFC& ufc, const Cell& cell, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains)
    
        


    .. cpp:function:: static void assemble_cell(arma::mat& A, UFC& ufc, const Cell& cell, const MeshFunction<uint>* domains)
    
        


    .. cpp:function:: static void assemble_exterior_facet(arma::mat& A, UFC& ufc, const Cell& cell, const Facet& facet, const uint local_facet, const MeshFunction<uint>* domains)
    
        


    .. cpp:function:: static void assemble_interior_facet(arma::mat& A, UFC& ufc, const Cell& cell, const Facet& facet, const uint local_facet, const MeshFunction<uint>* domains)
    
        


