
.. Documentation for the header file dolfin/fem/LocalSolver.h

.. _programmers_reference_cpp_fem_localsolver:

LocalSolver.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalSolver

    .. cpp:function:: void solve(GenericVector& x, const Form& a, const Form& L, bool symmetric=false) const
    
        Solve local (cell-wise) problem and copy result into global
        vector x.


