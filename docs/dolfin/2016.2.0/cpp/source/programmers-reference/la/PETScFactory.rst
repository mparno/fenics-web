
.. Documentation for the header file dolfin/la/PETScFactory.h

.. _programmers_reference_cpp_la_petscfactory:

PETScFactory.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PETScFactory

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearAlgebraFactory`
        
    .. cpp:function:: std::shared_ptr<GenericMatrix> create_matrix(MPI_Comm comm) const
    
        Create empty matrix


    .. cpp:function:: std::shared_ptr<GenericVector> create_vector(MPI_Comm comm) const
    
        Create empty vector


    .. cpp:function:: std::shared_ptr<TensorLayout> create_layout(std::size_t rank) const
    
        Create empty tensor layout


    .. cpp:function:: std::shared_ptr<GenericLinearOperator> create_linear_operator(MPI_Comm comm) const
    
        Create empty linear operator


    .. cpp:function:: std::shared_ptr<GenericLUSolver> create_lu_solver(MPI_Comm comm, std::string method) const
    
        Create LU solver


    .. cpp:function:: std::shared_ptr<GenericLinearSolver> create_krylov_solver(MPI_Comm comm, std::string method, std::string preconditioner) const
    
        Create Krylov solver


    .. cpp:function:: std::map<std::string, std::string> lu_solver_methods() const
    
        Return a list of available LU solver methods


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_methods() const
    
        Return a list of available Krylov solver methods


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_preconditioners() const
    
        Return a list of available preconditioners


    .. cpp:function:: static PETScFactory& instance()
    
        Return singleton instance


    .. cpp:function:: PETScFactory()
    
        Private constructor


