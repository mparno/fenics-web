
.. Documentation for the header file dolfin/la/uBLASILUPreconditioner.h

.. _programmers_reference_cpp_la_ublasilupreconditioner:

uBLASILUPreconditioner.h
========================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: uBLASILUPreconditioner

    *Parent class(es)*
    
        * :cpp:class:`uBLASPreconditioner`
        
    This class implements an incomplete LU factorization (ILU)
    preconditioner for the uBLAS Krylov solver.


    .. cpp:function:: uBLASILUPreconditioner(const Parameters& krylov_parameters)
    
        Constructor


    .. cpp:function:: void solve(uBLASVector& x, const uBLASVector& b) const
    
        Solve linear system Ax = b approximately


