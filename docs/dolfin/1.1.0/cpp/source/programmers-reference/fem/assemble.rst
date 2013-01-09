
.. Documentation for the header file dolfin/fem/assemble.h

.. _programmers_reference_cpp_fem_assemble:

assemble.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void assemble(GenericTensor& A, const Form& a, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble tensor


.. cpp:function:: void assemble(GenericTensor& A, const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble tensor on sub domain


.. cpp:function:: void assemble(GenericTensor& A, const Form& a, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble tensor on sub domains


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble system (A, b)


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const DirichletBC& bc, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble system (A, b) and apply Dirichlet boundary condition


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*> bcs, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble system (A, b) and apply Dirichlet boundary conditions


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*> bcs, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains, const GenericVector* x0, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Assemble system (A, b) on sub domains and apply Dirichlet boundary conditions


.. cpp:function:: void symmetric_assemble(GenericMatrix& As, GenericMatrix& An, const Form& a, const std::vector<const DirichletBC*> bcs, const MeshFunction<std::size_t>* cell_domains=NULL, const MeshFunction<std::size_t>* exterior_facet_domains=NULL, const MeshFunction<std::size_t>* interior_facet_domains=NULL, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Symmetric assembly of As, storing the modifications in An. To create
    matching RHS, assemble and apply bcs normally, then subtract An*b.
    In this variant of symmetric_assemble, rows and columns use the same BCs.


.. cpp:function:: void symmetric_assemble(GenericMatrix& As, GenericMatrix& An, const Form& a, const std::vector<const DirichletBC*> row_bcs, const std::vector<const DirichletBC*> col_bcs, const MeshFunction<std::size_t>* cell_domains=NULL, const MeshFunction<std::size_t>* exterior_facet_domains=NULL, const MeshFunction<std::size_t>* interior_facet_domains=NULL, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true, bool keep_diagonal=false)

    Symmetric assembly of As, storing the modifications in An. To create
    matching RHS, assemble and apply bcs normally, then subtract An*b.
    In this variant of symmetric_assemble, rows and columns use (potentially)
    different BCs. The BCs will be different for example in coupling
    (i.e., off-diagonal) blocks of a block matrix.


.. cpp:function:: double assemble(const Form& a, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)

    Assemble scalar


.. cpp:function:: double assemble(const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)

    Assemble scalar on sub domain


.. cpp:function:: double assemble(const Form& a, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)

    Assemble scalar on sub domains


