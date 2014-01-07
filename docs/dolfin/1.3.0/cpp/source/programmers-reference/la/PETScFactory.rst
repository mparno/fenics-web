
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
        
    .. cpp:function:: boost::shared_ptr<GenericMatrix> create_matrix() const
    
        Create empty matrix


    .. cpp:function:: boost::shared_ptr<GenericVector> create_vector() const
    
        Create empty vector (global)


    .. cpp:function:: boost::shared_ptr<GenericVector> create_local_vector() const
    
        Create empty vector (local)


    .. cpp:function:: boost::shared_ptr<TensorLayout> create_layout(std::size_t rank) const
    
        Create empty tensor layout


    .. cpp:function:: boost::shared_ptr<GenericLinearOperator> create_linear_operator() const
    
        Create empty linear operator


    .. cpp:function:: boost::shared_ptr<GenericLUSolver> create_lu_solver(std::string method) const
    
        Create LU solver


    .. cpp:function:: boost::shared_ptr<GenericLinearSolver> create_krylov_solver(std::string method, std::string preconditioner) const
    
        Create Krylov solver


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > lu_solver_methods() const
    
        Return a list of available LU solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_methods() const
    
        Return a list of available Krylov solver methods


    .. cpp:function:: std::vector<std::pair<std::string, std::string> > krylov_solver_preconditioners() const
    
        Return a list of available preconditioners


    .. cpp:function:: static PETScFactory& instance()
    
        Return singleton instance


    .. cpp:function:: PETScFactory()
    
        Private constructor


