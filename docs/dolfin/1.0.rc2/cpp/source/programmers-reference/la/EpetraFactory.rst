
.. Documentation for the header file dolfin/la/EpetraFactory.h

.. _programmers_reference_cpp_la_epetrafactory:

EpetraFactory.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EpetraFactory

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: EpetraMatrix* create_matrix() const
    
        --- LinearAlgebraFactory interface
        Create empty matrix


    .. cpp:function:: EpetraVector* create_vector() const
    
        Create empty vector (global)


    .. cpp:function:: EpetraVector* create_local_vector() const
    
        Create empty vector (local)


    .. cpp:function:: SparsityPattern* create_pattern() const
    
        Create empty sparsity pattern


    .. cpp:function:: EpetraLUSolver* create_lu_solver(std::string method) const
    
        Create LU solver


    .. cpp:function:: EpetraKrylovSolver* create_krylov_solver(std::string method, std::string preconditioner) const
    
        Create Krylov solver


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > lu_solver_methods() const
    
        Return a list of available LU solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_methods() const
    
        Return a list of available Krylov solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_preconditioners() const
    
        Return a list of available preconditioners


