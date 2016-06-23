
.. Documentation for the header file dolfin/la/Ifpack2Preconditioner.h

.. _programmers_reference_cpp_la_ifpack2preconditioner:

Ifpack2Preconditioner.h
=======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Ifpack2Preconditioner

    *Parent class(es)*
    
        * :cpp:class:`TrilinosPreconditioner`
        
        * :cpp:class:`Variable`
        
    Implements preconditioners using Ifpack2 from Trilinos


    .. cpp:function:: explicit Ifpack2Preconditioner(std::string type = "default")
    
        Create a particular preconditioner object


    .. cpp:function:: void set(BelosKrylovSolver& solver)
    
        Set the preconditioner type on a solver


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


    .. cpp:function:: void init(std::shared_ptr<const TpetraMatrix> P)
    
        Initialise preconditioner based on Operator P


    .. cpp:function:: static std::map<std::string, std::string> preconditioners()
    
        Return a list of available preconditioners


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


