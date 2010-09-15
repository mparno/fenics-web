.. Documentation for the header file dolfin/la/PETScPreconditioner.h

.. _programmers_reference_cpp_la_petscpreconditioner:

PETScPreconditioner.h
=====================

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: PETScPreconditioner

    *Parent class*
    
        * :cpp:class:`PETScObject,`
        
    This class is a wrapper for configuring PETSc preconditioners. It does
    not own a preconditioner. It can take a PETScKrylovSolver and set the
    preconditioner type and parameters.

    .. cpp:function:: explicit PETScPreconditioner(std::string type = "default")
    
        Create a particular preconditioner object

    .. cpp:function:: void set(PETScKrylovSolver& solver) const
    
        Set the precondtioner type and parameters

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values

    .. cpp:function:: std::string type
    
        Named preconditioner

