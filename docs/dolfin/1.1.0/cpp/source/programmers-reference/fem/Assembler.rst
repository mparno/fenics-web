
.. Documentation for the header file dolfin/fem/Assembler.h

.. _programmers_reference_cpp_fem_assembler:

Assembler.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Assembler

    *Parent class(es)*
    
        * :cpp:class:`AssemblerBase`
        
    This class provides automated assembly of linear systems, or
    more generally, assembly of a sparse tensor from a given
    variational form.
    
    Subdomains for cells and facets may be specified in a number of
    different ways:
    
    1. By explicitly passing :cpp:class:`MeshFunction` (as pointers) to the
       assemble functions
    
    2. By assigning subdomain indicators specified by :cpp:class:`MeshFunction`
       to the :cpp:class:`Form` being assembled:
    
       .. code-block:: c++
    
           form.dx = cell_domains
           form.ds = exterior_facet_domains
           form.dS = interior_facet_domains
    
    3. By markers stored as part of the :cpp:class:`Mesh` (in :cpp:class:`MeshDomains`)
    
    4. By specifying a :cpp:class:`SubDomain` which specifies the domain numbered
       as 0 (with the rest treated as domain number 1)
    
    Note that (1) overrides (2), which overrides (3).


    .. cpp:function:: void assemble(GenericTensor& A, const Form& a)
    
        Assemble tensor from given form
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.


    .. cpp:function:: void assemble(GenericTensor& A, const Form& a, const SubDomain& sub_domain)
    
        Assemble tensor from given form on subdomain
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain to assemble on.


    .. cpp:function:: void assemble(GenericTensor& A, const Form& a, const MeshFunction<std::size_t>* cell_domains, const MeshFunction<std::size_t>* exterior_facet_domains, const MeshFunction<std::size_t>* interior_facet_domains)
    
        Assemble tensor from given form on subdomains
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.
            cell_domains (:cpp:class:`MeshFunction` <std::size_t>)
                Cell domains.
            exterior_facet_domains (:cpp:class:`MeshFunction` <std::size_t>)
                The exterior facet domains.
            interior_facet_domains (:cpp:class:`MeshFunction` <std::size_t>)
                The interior facet domains.


    .. cpp:function:: void assemble_cells(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<std::size_t>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over cells. This function is
        provided for users who wish to build a customized assembler.


    .. cpp:function:: void assemble_exterior_facets(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<std::size_t>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over exterior facets. This
        function is provided for users who wish to build a customized
        assembler.


    .. cpp:function:: void assemble_interior_facets(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<std::size_t>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over interior facets. This
        function is provided for users who wish to build a customized
        assembler.


    .. cpp:function:: void add_to_global_tensor(GenericTensor& A, std::vector<double>& cell_tensor, std::vector<const std::vector<dolfin::la_index>* >& dofs)
    
        Add cell tensor to global tensor. Hook to allow the SymmetricAssembler
        to split the cell tensor into symmetric/antisymmetric parts.


