.. Documentation for the header file dolfin/fem/assemble.h

.. _programmers_reference_cpp_fem_assemble:

assemble.h
==========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

    .. cpp:function:: double assemble(const Form& a, bool reset_sparsity=true, bool add_values=false)
    
        Assemble scalar

    .. cpp:function:: double assemble(const Form& a, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false)
    
        Assemble scalar on sub domains

    .. cpp:function:: double assemble(const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false)
    
        Assemble scalar on sub domain

    .. cpp:function:: void assemble(GenericTensor& A, const Form& a, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor

    .. cpp:function:: void assemble(GenericTensor& A, const Form& a, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor on sub domains

    .. cpp:function:: void assemble(GenericTensor& A, const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor on sub domain

    .. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, bool reset_sparsitys=true, bool add_values=false)
    
        Assemble system (A, b)

    .. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const DirichletBC& bc, bool reset_sparsitys=true, bool add_values=false)
    
        Assemble system (A, b) and apply Dirichlet boundary condition

    .. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*>& bcs, bool reset_sparsitys=true, bool add_values=false)
    
        Assemble system (A, b) and apply Dirichlet boundary conditions

    .. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*>& bcs, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, const GenericVector* x0, bool reset_sparsitys=true, bool add_values=false)
    
        Assemble system (A, b) on sub domains and apply Dirichlet boundary conditions

