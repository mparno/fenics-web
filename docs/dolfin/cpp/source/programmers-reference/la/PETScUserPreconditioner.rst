
.. Documentation for the header file dolfin/la/PETScUserPreconditioner.h

.. _programmers_reference_cpp_la_petscuserpreconditioner:

PETScUserPreconditioner.h
=========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScUserPreconditioner

    *Parent class(es)*
    
        * :cpp:class:`PETScObject`
        
    This class specifies the interface for user-defined Krylov
    method PETScPreconditioners. A user wishing to implement her own
    PETScPreconditioner needs only supply a function that approximately
    solves the linear system given a right-hand side.


    .. cpp:function:: PETScUserPreconditioner()
    
        Constructor


    .. cpp:function:: void solve(PETScVector& x, const PETScVector& b) = 0
    
        Solve linear system approximately for given right-hand side b


