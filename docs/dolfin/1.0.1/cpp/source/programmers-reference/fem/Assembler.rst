
.. Documentation for the header file dolfin/fem/Assembler.h

.. _programmers_reference_cpp_fem_assembler:

Assembler.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Assembler

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


    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)
    
        Assemble tensor from given form
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.
            reset_sparsity (bool)
                Optional argument: Default value is true.
                This controls whether the sparsity pattern of the
                given tensor is reset prior to assembly.
            add_values (bool)
                Optional argument: Default value is false.
                This controls whether values are added to the given
                tensor or if it is zeroed prior to assembly.
            finalize_tensor (bool)
                Optional argument: Default value is true.
                This controls whether the assembler finalizes the
                given tensor after assembly is completed by calling
                A.apply().


    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)
    
        Assemble tensor from given form on subdomain
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain to assemble on.
            reset_sparsity (bool)
                Optional argument: Default value is true.
                This controls whether the sparsity pattern of the
                given tensor is reset prior to assembly.
            add_values (bool)
                Optional argument: Default value is false.
                This controls whether values are added to the given
                tensor or if it is zeroed prior to assembly.
            finalize_tensor (bool)
                Optional argument: Default value is true.
                This controls whether the assembler finalizes the
                given tensor after assembly is completed by calling
                A.apply().


    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false, bool finalize_tensor=true)
    
        Assemble tensor from given form on subdomains
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.
            cell_domains (:cpp:class:`MeshFunction` <uint>)
                Cell domains.
            exterior_facet_domains (:cpp:class:`MeshFunction` <uint>)
                The exterior facet domains.
            interior_facet_domains (:cpp:class:`MeshFunction` <uint>)
                The interior facet domains.
            reset_sparsity (bool)
                Optional argument: Default value is true.
                This controls whether the sparsity pattern of the
                given tensor is reset prior to assembly.
            add_values (bool)
                Optional argument: Default value is false.
                This controls whether values are added to the given
                tensor or if it is zeroed prior to assembly.
            finalize_tensor (bool)
                Optional argument: Default value is true.
                This controls whether the assembler finalizes the
                given tensor after assembly is completed by calling
                A.apply().


    .. cpp:function:: static void assemble_cells(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<uint>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over cells. This function is
        provided for users who wish to build a customized assembler.


    .. cpp:function:: static void assemble_exterior_facets(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<uint>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over exterior facets. This
        function is provided for users who wish to build a customized
        assembler.


    .. cpp:function:: static void assemble_interior_facets(GenericTensor& A, const Form& a, UFC& ufc, const MeshFunction<uint>* domains, std::vector<double>* values)
    
        Assemble tensor from given form over interior facets. This
        function is provided for users who wish to build a customized
        assembler.


