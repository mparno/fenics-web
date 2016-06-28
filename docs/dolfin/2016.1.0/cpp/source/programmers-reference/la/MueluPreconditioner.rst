
.. Documentation for the header file dolfin/la/MueluPreconditioner.h

.. _programmers_reference_cpp_la_muelupreconditioner:

MueluPreconditioner.h
=====================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MueluPreconditioner

    *Parent class(es)*
    
        * :cpp:class:`TrilinosPreconditioner`
        
        * :cpp:class:`Variable`
        
    Implements Muelu preconditioner from Trilinos


    .. cpp:function:: MueluPreconditioner()
    
        Create a particular preconditioner object


    .. cpp:function:: void set(BelosKrylovSolver& solver)
    
        Set the preconditioner on a solver


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void init(std::shared_ptr<const TpetraMatrix> P)
    
        Initialise preconditioner based on Operator P


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


