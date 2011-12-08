
.. Documentation for the header file dolfin/la/MTL4Factory.h

.. _programmers_reference_cpp_la_mtl4factory:

MTL4Factory.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MTL4Factory

    *Parent class(es)*
    
        * :cpp:class:`LinearAlgebraFactory`
        
    .. cpp:function:: MTL4Matrix* create_matrix() const
    
        Create empty matrix


    .. cpp:function:: MTL4Vector* create_vector() const
    
        Create empty vector (global)


    .. cpp:function:: MTL4Vector* create_local_vector() const
    
        Create empty vector (local)


    .. cpp:function:: GenericSparsityPattern* create_pattern() const
    
        Dummy sparsity pattern


    .. cpp:function:: UmfpackLUSolver* create_lu_solver(std::string method) const
    
        Create LU solver


    .. cpp:function:: ITLKrylovSolver* create_krylov_solver(std::string method, std::string preconditioner) const
    
        Create Krylov solver


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > lu_solver_methods() const
    
        Return a list of available LU solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_methods() const
    
        Return a list of available Krylov solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_preconditioners() const
    
        Return a list of available preconditioners


