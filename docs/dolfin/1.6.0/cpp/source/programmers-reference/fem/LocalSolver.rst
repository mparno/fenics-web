
.. Documentation for the header file dolfin/fem/LocalSolver.h

.. _programmers_reference_cpp_fem_localsolver:

LocalSolver.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: LocalSolver

    .. cpp:function:: LocalSolver(std::shared_ptr<const Form> a, std::shared_ptr<const Form> L, SolverType solver_type=LU)
    
        Constructor (shared pointer version)


    .. cpp:function:: LocalSolver(std::shared_ptr<const Form> a, SolverType solver_type=LU)
    
        Constructor (shared pointer version)


    .. cpp:function:: void solve_global_rhs(Function& u) const
    
        Solve local (cell-wise) problems A_e x_e = b_e, where A_e is
        the cell matrix LHS and b_e is the global RHS vector b
        restricted to the cell, i.e. b_e may contain contributions
        from neighbouring cells. The solution is exact for the case in
        which there is no coupling between cell contributions to the
        global matrix A, e.g. the discontinuous Galerkin matrix. The
        result is copied into x.


    .. cpp:function:: void solve_local_rhs(Function& u) const
    
        Solve local (cell-wise) problems A_e x_e = b_e where A_e and
        b_e are the cell element tensors. This function is useful for
        computing (approximate) cell-wise projections, for example for
        post-processing. It much more efficient than computing global
        projections.


    .. cpp:function:: void solve_local(GenericVector& x, const GenericVector& b, const GenericDofMap& dofmap_b) const
    
        Solve local problems for given RHS and corresponding dofmap
        for RHS


    .. cpp:function:: void factorize()
    
        Factorise LHS for all cells and store


    .. cpp:function:: void clear_factorization()
    
        Reset (clear) any stored factorisations


