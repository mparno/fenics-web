
.. Documentation for the header file dolfin/fem/CCFEMAssembler.h

.. _programmers_reference_cpp_fem_ccfemassembler:

CCFEMAssembler.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: CCFEMAssembler

    This class implements functionality for finite element assembly
    over cut and composite finite element (CCFEM) function spaces.


    .. cpp:function:: CCFEMAssembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericTensor& A, const CCFEMForm& a)
    
        Assemble tensor from given form
        
        *Arguments*
            A (:cpp:class:`GenericTensor`)
                The tensor to assemble.
            a (:cpp:class:`Form`)
                The form to assemble the tensor from.


