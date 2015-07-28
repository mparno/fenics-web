
.. Documentation for the header file dolfin/adaptivity/LocalAssembler.h

.. _programmers_reference_cpp_adaptivity_localassembler:

LocalAssembler.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalAssembler

    


    .. cpp:function:: static void assemble(Eigen::MatrixXd& A, UFC& ufc, const std::vector<double>& vertex_coordinates, ufc::cell& ufc_cell, const Cell& cell, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains)
    
        


    .. cpp:function:: static void assemble_cell(Eigen::MatrixXd& A, UFC& ufc, const std::vector<double>& vertex_coordinates, const ufc::cell& ufc_cell, const Cell& cell, const MeshFunction<std::size_t>* domains)
    
        


    .. cpp:function:: static void assemble_exterior_facet(Eigen::MatrixXd& A, UFC& ufc, const std::vector<double>& vertex_coordinates, const ufc::cell& ufc_cell, const Cell& cell, const Facet& facet, const std::size_t local_facet, const MeshFunction<std::size_t>* domains)
    
        


    .. cpp:function:: static void assemble_interior_facet(Eigen::MatrixXd& A, UFC& ufc, const std::vector<double>& vertex_coordinates, const ufc::cell& ufc_cell, const Cell& cell, const Facet& facet, const std::size_t local_facet, const MeshFunction<std::size_t>* domains)
    
        


