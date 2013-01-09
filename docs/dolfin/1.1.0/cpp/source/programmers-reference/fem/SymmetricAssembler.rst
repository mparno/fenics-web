
.. Documentation for the header file dolfin/fem/SymmetricAssembler.h

.. _programmers_reference_cpp_fem_symmetricassembler:

SymmetricAssembler.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SymmetricAssembler

    *Parent class(es)*
    
        * :cpp:class:`Assembler`
        
    This class provides implements an assembler for systems
    of the form Ax = b. Its assembly algorithms are similar to SystemAssember's,
    but it saves the matrix modifications into a separate tensor so that it
    can later apply the symmetric modifications to a RHS vector.
    The non-symmetric part is only nonzero in BC columns, and is zero in all BC
    rows, so that [(A_s+A_n) x = b] implies [A_s x = b - A_n b], IF b has
    boundary conditions applied. (If the final A is composed from a sum of
    A_s matrices, BCs must be reapplied to make the diagonal value for BC
    dofs 1.0. The matrix will remain symmetric since only the diagonal is
    changed.)
    
    *Example*
    
       .. code-block:: c++
    
           std::vector<const DirichletBC*> bcs = {bc};
           SymmetricAssembler::assemble(A, A_n, a, bcs, bcs);
           Assembler::assemble(b, L);
           bc.apply(b)
           A_n.mult(b, b_mod);
           b -= b_mod;


    .. cpp:function:: SymmetricAssembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericMatrix& A, GenericMatrix& B, const Form& a, const std::vector<const DirichletBC*> row_bcs, const std::vector<const DirichletBC*> col_bcs, const MeshFunction<std::size_t>* cell_domains=NULL, const MeshFunction<std::size_t>* exterior_facet_domains=NULL, const MeshFunction<std::size_t>* interior_facet_domains=NULL)
    
        Assemble a and apply Dirichlet boundary conditions. Returns two
        matrices, where the second contains the symmetric modifications
        suitable for modifying RHS vectors.
        
        Note: row_bcs and col_bcs will normally be the same, but are different
        for e.g. off-diagonal block matrices in a mixed PDE.


    .. cpp:function:: void assemble(GenericMatrix& A, GenericMatrix& B, const Form& a, const std::vector<const DirichletBC*> row_bcs, const std::vector<const DirichletBC*> col_bcs, const SubDomain& sub_domain)
    
        Assemble a and apply Dirichlet boundary conditions. Returns two
        matrices, where the second contains the symmetric modifications
        suitable for modifying RHS vectors.
        
        Note: row_bcs and col_bcs will normally be the same, but are different
        for e.g. off-diagonal block matrices in a mixed PDE.


    .. cpp:function:: void add_to_global_tensor(GenericTensor& A, std::vector<double>& local_A, std::vector<const std::vector<dolfin::la_index>* >& dofs)
    
        Add cell tensor to global tensor. Hook to allow the SymmetricAssembler
        to split the cell tensor into symmetric/antisymmetric parts.


