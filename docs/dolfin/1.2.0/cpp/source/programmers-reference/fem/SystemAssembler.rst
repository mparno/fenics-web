
.. Documentation for the header file dolfin/fem/SystemAssembler.h

.. _programmers_reference_cpp_fem_systemassembler:

SystemAssembler.h
=================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SystemAssembler

    *Parent class(es)*
    
        * :cpp:class:`AssemblerBase`
        
    This class provides implements an assembler for systems
    of the form Ax = b. It differs from the default DOLFIN
    assembler in that it assembles both A and b and the same
    time (leading to better performance) and in that it applies
    boundary conditions at the time of assembly.


    .. cpp:function:: SystemAssembler()
    
        Constructor


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L)
    
        Assemble system (A, b)


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const DirichletBC& bc)
    
        Assemble system (A, b) and apply Dirichlet boundary condition


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*> bcs)
    
        Assemble system (A, b) and apply Dirichlet boundary conditions


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, const std::vector<const DirichletBC*> bcs, const GenericVector* x0)
    
        Assemble system (A, b) and apply Dirichlet boundary conditions


.. cpp:class:: Scratch

