
.. Documentation for the header file dolfin/fem/LocalAssembler.h

.. _programmers_reference_cpp_fem_localassembler:

LocalAssembler.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalAssembler

    Assembly of local cell tensors. Used by the adaptivity and
    LocalSolver functionality in dolfin. The local assembly
    functionality provided here is also wrapped as a free function
    assemble_local(form_a, cell) in Python for easier usage. Used
    from C++ the versions defined below will be faster than the free
    function as less objects need to be created and thrown away


    .. cpp:function:: static void assemble(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A, UFC& ufc, const std::vector<double>& coordinate_dofs, ufc::cell& ufc_cell, const Cell& cell, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains)
    
        Assemble a local tensor on a cell


    .. cpp:function:: static void assemble_cell(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A, UFC& ufc, const std::vector<double>& coordinate_dofs, const ufc::cell& ufc_cell, const Cell& cell, const MeshFunction<std::size_t>* cell_domains)
    
        Worker method called by assemble() to perform assembly of
        volume integrals (dx)


    .. cpp:function:: static void assemble_exterior_facet(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A, UFC& ufc, const std::vector<double>& coordinate_dofs, const ufc::cell& ufc_cell, const Cell& cell, const Facet& facet, const std::size_t local_facet, const MeshFunction<std::size_t>* exterior_facet_domains)
    
        Worker method called by assemble() for each of the cell's
        external facets to perform assembly of external facet
        integrals (ds)


    .. cpp:function:: static void assemble_interior_facet(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& A, UFC& ufc, const std::vector<double>& coordinate_dofs, const ufc::cell& ufc_cell, const Cell& cell, const Facet& facet, const std::size_t local_facet, const MeshFunction<std::size_t>* interior_facet_domains, const MeshFunction<std::size_t>* cell_domains)
    
        Worker method called by assemble() for each of the cell's
        internal facets to perform assembly of internal facet
        integrals (dS)


