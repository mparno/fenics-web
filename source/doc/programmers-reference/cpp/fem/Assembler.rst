.. Documentation for the header file dolfin/fem/Assembler.h

.. _programmers_reference_cpp_fem_assembler:

Assembler.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Assembler

    This class provides automated assembly of linear systems, or
    more generally, assembly of a sparse tensor from a given
    variational form.
    
    Subdomains for cells and facets may be specified in a number
    of different ways:
    
    1. By explicitly passing MeshFunctions (as pointers) to the
       assemble functions
    
    2. By assigning subdomain indicators specified by MeshFunctions
       to the Form being assembled:
    
       form.cell_domains = cell_domains
       form.exterior_facet_domains = exterior_facet_domains
       form.interior_facet_domains = interior_facet_domains
    
    3. By MeshFunctions stored in MeshData as
    
       "cell_domains"
       "exterior_facet_domains"
       "interior_facet_domains"
    
    4. By specifying a SubDomain which specifies the domain numbered
       as 0 (with the rest treated as domain number 1)
    
    Note that (1) overrides (2), which overrides (3).

    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor from given form

    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, const SubDomain& sub_domain, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor from given form on sub domain

    .. cpp:function:: static void assemble(GenericTensor& A, const Form& a, const MeshFunction<uint>* cell_domains, const MeshFunction<uint>* exterior_facet_domains, const MeshFunction<uint>* interior_facet_domains, bool reset_sparsity=true, bool add_values=false)
    
        Assemble tensor from given form on sub domains

