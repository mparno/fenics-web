
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
    
    Subdomains for cells and facets may be specified by assigning
    subdomain indicators specified by :cpp:class:`MeshFunction` to the :cpp:class:`Form`
    being assembled:
    
       .. code-block:: c++
    
           form.dx = cell_domains
           form.ds = exterior_facet_domains
           form.dS = interior_facet_domains


    .. cpp:function:: Assembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericTensor& A, const Form& a)
    
        Assemble tensor from given form
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.


    .. cpp:function:: void assemble_cells(GenericTensor& A, const Form& a, UFC& ufc, std::shared_ptr<const MeshFunction<std::size_t>> domains, std::vector<double>* values)
    
        Assemble tensor from given form over cells. This function is
        provided for users who wish to build a customized assembler.


    .. cpp:function:: void assemble_exterior_facets(GenericTensor& A, const Form& a, UFC& ufc, std::shared_ptr<const MeshFunction<std::size_t>> domains, std::vector<double>* values)
    
        Assemble tensor from given form over exterior facets. This
        function is provided for users who wish to build a customized
        assembler.


    .. cpp:function:: void assemble_interior_facets(GenericTensor& A, const Form& a, UFC& ufc, std::shared_ptr<const MeshFunction<std::size_t>> domains, std::shared_ptr<const MeshFunction<std::size_t>> cell_domains, std::vector<double>* values)
    
        Assemble tensor from given form over interior facets. This
        function is provided for users who wish to build a customized
        assembler.


    .. cpp:function:: void assemble_vertices(GenericTensor& A, const Form& a, UFC& ufc, std::shared_ptr<const MeshFunction<std::size_t>> domains)
    
        Assemble tensor from given form over vertices. This function is
        provided for users who wish to build a customized assembler.


