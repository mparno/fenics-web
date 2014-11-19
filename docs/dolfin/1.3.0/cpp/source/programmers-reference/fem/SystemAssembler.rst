
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
        
    This class provides an assembler for systems of the form
    Ax = b. It differs from the default DOLFIN assembler in that it
    applies boundary conditions at the time of assembly, which
    preserves any symmetries in A.


    .. cpp:function:: SystemAssembler(const Form& a, const Form& L)
    
        Constructor


    .. cpp:function:: SystemAssembler(const Form& a, const Form& L, const DirichletBC& bc)
    
        Constructor


    .. cpp:function:: SystemAssembler(const Form& a, const Form& L, const std::vector<const DirichletBC*> bcs)
    
        Constructor


    .. cpp:function:: SystemAssembler(boost::shared_ptr<const Form> a, boost::shared_ptr<const Form> L)
    
        Constructor


    .. cpp:function:: SystemAssembler(boost::shared_ptr<const Form> a, boost::shared_ptr<const Form> L, const DirichletBC& bc)
    
        Constructor


    .. cpp:function:: SystemAssembler(boost::shared_ptr<const Form> a, boost::shared_ptr<const Form> L, const std::vector<const DirichletBC*> bcs)
    
        Constructor


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b)
    
        Assemble system (A, b)


    .. cpp:function:: void assemble(GenericMatrix& A)
    
        Assemble matrix A


    .. cpp:function:: void assemble(GenericVector& b)
    
        Assemble vector b


    .. cpp:function:: void assemble(GenericMatrix& A, GenericVector& b, const GenericVector& x0)
    
        Assemble system (A, b) for (negative) increment dx, where
        x = x0 - dx is solution to system a == -L subject to bcs.
        Suitable for use inside a (quasi-)Newton solver.


    .. cpp:function:: void assemble(GenericVector& b, const GenericVector& x0)
    
        Assemble rhs vector b for (negative) increment dx, where
        x = x0 - dx is solution to system a == -L subject to bcs.
        Suitable for use inside a (quasi-)Newton solver.


.. cpp:class:: Scratch

