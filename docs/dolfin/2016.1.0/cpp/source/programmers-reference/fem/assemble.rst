
.. Documentation for the header file dolfin/fem/assemble.h

.. _programmers_reference_cpp_fem_assemble:

assemble.h
==========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    



.. cpp:function:: void assemble(GenericTensor& A, const Form& a)

    Assemble tensor


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, std::vector<std::shared_ptr<const DirichletBC>> bcs)

    Assemble system (A, b) and apply Dirichlet boundary conditions


.. cpp:function:: void assemble_system(GenericMatrix& A, GenericVector& b, const Form& a, const Form& L, std::vector<std::shared_ptr<const DirichletBC>> bcs, const GenericVector& x0)

    Assemble system (A, b) on sub domains and apply Dirichlet
    boundary conditions


.. cpp:function:: void assemble_multimesh(GenericTensor& A, const MultiMeshForm& a)

    Assemble tensor from multimesh form


.. cpp:function:: double assemble(const Form& a)

    Assemble scalar


