
.. Documentation for the header file dolfin/fem/OpenMpAssembler.h

.. _programmers_reference_cpp_fem_openmpassembler:

OpenMpAssembler.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: OpenMpAssembler

    *Parent class(es)*
    
        * :cpp:class:`AssemblerBase`
        
    This class provides automated assembly of linear systems, or
    more generally, assembly of a sparse tensor from a given
    variational form.
    
    The MeshFunction arguments can be used to specify assembly over
    subdomains of the mesh cells, exterior facets or interior
    facets. Either a null pointer or an empty MeshFunction may be
    used to specify that the tensor should be assembled over the
    entire set of cells or facets.


    .. cpp:function:: OpenMpAssembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericTensor& A, const Form& a)
    
        Assemble tensor from given form


