
.. Documentation for the header file dolfin/fem/SystemAssembler.h

.. _programmers_reference_cpp_fem_systemassembler:

SystemAssembler.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SystemAssembler

    This class provides implements an assembler for systems
    of the form Ax = b. It differs from the default DOLFIN
    assembler in that it assembles both A and b and the same
    time (leading to better performance) and in that it applies
    boundary conditions at the time of assembly.


    .. cpp:function:: static void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, bool reset_sparsity=true, bool add_values=false)
    
        Assemble system (A, b)


    .. cpp:function:: static void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const DirichletBC& bc, bool reset_sparsity=true, bool add_values=true)
    
        Assemble system (A, b) and apply Dirichlet boundary condition


    .. cpp:function:: static void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*>& bcs, bool reset_sparsity=true, bool add_values=false)
    
        Assemble system (A, b) and apply Dirichlet boundary conditions


    .. cpp:function:: static void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*>& bcs, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, const GenericVector* x0, bool reset_sparsity=true, bool add_values=false)
    
        Assemble system (A, b) and apply Dirichlet boundary conditions


.. cpp:class:: Scratch

