
.. Documentation for the header file dolfin/la/solve.h

.. _programmers_reference_cpp_la_solve:

solve.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: std::size_t solve(const GenericLinearOperator& A, GenericVector& x, const GenericVector& b, std::string method = "lu", std::string preconditioner = "none")

    Solve linear system Ax = b


.. cpp:function:: void list_linear_algebra_backends()

    List available linear algebra backends


.. cpp:function:: void list_linear_solver_methods()

    List available solver methods for current linear algebra backend


.. cpp:function:: void list_lu_solver_methods()

    List available LU methods for current linear algebra backend


.. cpp:function:: void list_krylov_solver_methods()

    List available Krylov methods for current linear algebra backend


.. cpp:function:: void list_krylov_solver_preconditioners()

    List available preconditioners for current linear algebra
    backend


.. cpp:function:: bool has_lu_solver_method(std::string method)

    Return true if LU method for the current linear algebra backend is
    available


.. cpp:function:: bool has_krylov_solver_method(std::string method)

    Return true if Krylov method for the current linear algebra
    backend is available


.. cpp:function:: bool has_krylov_solver_preconditioner(std::string preconditioner)

    Return true if Preconditioner for the current linear algebra
    backend is available


.. cpp:function:: std::vector<std::pair<std::string, std::string> > linear_algebra_backends()

    Return available linear algebra backends


.. cpp:function:: std::vector<std::pair<std::string, std::string> > linear_solver_methods()

    Return a list of available solver methods for current linear
    algebra backend


.. cpp:function:: std::vector<std::pair<std::string, std::string> > lu_solver_methods()

    Return a list of available LU methods for current linear algebra
    backend


.. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_methods()

    Return a list of available Krylov methods for current linear
    algebra backend


.. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_preconditioners()

    Return a list of available preconditioners for current linear
    algebra backend


.. cpp:function:: double residual(const GenericLinearOperator& A, const GenericVector& x, const GenericVector& b)

    Compute residual ||Ax - b||


.. cpp:function:: double norm(const GenericVector& x, std::string norm_type="l2")

    Compute norm of vector. Valid norm types are "l2", "l1" and
    "linf".


.. cpp:function:: double normalize(GenericVector& x, std::string normalization_type = "average")

    Normalize vector according to given normalization type


