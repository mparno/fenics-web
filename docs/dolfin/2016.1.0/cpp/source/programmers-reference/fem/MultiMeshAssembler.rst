
.. Documentation for the header file dolfin/fem/MultiMeshAssembler.h

.. _programmers_reference_cpp_fem_multimeshassembler:

MultiMeshAssembler.h
====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshAssembler

    *Parent class(es)*
    
        * :cpp:class:`AssemblerBase`
        
    This class implements functionality for finite element assembly
    over cut and composite finite element (MultiMesh) function spaces.


    .. cpp:function:: MultiMeshAssembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericTensor& A, const MultiMeshForm& a)
    
        Assemble tensor from given form
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.


