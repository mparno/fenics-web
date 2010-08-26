.. Documentation for the header file dolfin/fem/Assembler.h

.. _programmers_reference_cpp_fem_Mesh:

Assembler.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Assembler

        This class provides automated assembly of linear systems, or
        more generally, assembly of a sparse tensor from a given
        variational form.
        
        The MeshFunction arguments can be used to specify assembly over
        subdomains of the mesh cells, exterior facets or interior
        facets. Either a null pointer or an empty MeshFunction may be
        used to specify that the tensor should be assembled over the
        entire set of cells or facets.

    .. cpp:function:: static void assemble(GenericTensor& A,
                                           const Form& a,
                                           const MeshFunction<uint>* cell_domains,
                                           const MeshFunction<uint>* exterior_facet_domains,
                                           const MeshFunction<uint>* interior_facet_domains,
                                           bool reset_sparsity=true,
                                           bool add_values=false)
    
        Assemble tensor on sub domains

    .. cpp:function:: static void assemble(GenericTensor& A,
                                           const Form& a,
                                           const SubDomain& sub_domain,
                                           bool reset_sparsity=true,
                                           bool add_values=false)
    
        Assemble tensor on sub domain

    .. cpp:function:: static void assemble(GenericTensor& A,
                       const Form& a,
                       bool reset_sparsity=true,
                       bool add_values=false)
    
        Assemble tensor

