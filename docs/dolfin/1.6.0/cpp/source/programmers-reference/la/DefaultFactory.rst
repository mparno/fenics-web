
.. Documentation for the header file dolfin/la/DefaultFactory.h

.. _programmers_reference_cpp_la_defaultfactory:

DefaultFactory.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DefaultFactory

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearAlgebraFactory`
        
    .. cpp:function:: DefaultFactory()
    
        Constructor


    .. cpp:function:: std::shared_ptr<GenericMatrix> create_matrix() const
    
        Create empty matrix


    .. cpp:function:: std::shared_ptr<GenericVector> create_vector() const
    
        Create empty vector


    .. cpp:function:: std::shared_ptr<TensorLayout> create_layout(std::size_t rank) const
    
        Create empty tensor layout


    .. cpp:function:: std::shared_ptr<GenericLinearOperator> create_linear_operator() const
    
        Create empty linear operator


    .. cpp:function:: std::shared_ptr<dolfin::GenericLUSolver> create_lu_solver(std::string method) const
    
        Create LU solver


    .. cpp:function:: std::shared_ptr<dolfin::GenericLinearSolver> create_krylov_solver(std::string method, std::string preconditioner) const
    
        Create Krylov solver


    .. cpp:function:: std::map<std::string, std::string> lu_solver_methods() const
    
        Return a list of available LU solver methods


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_methods() const
    
        Return a list of available Krylov solver methods


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_preconditioners() const
    
        Return a list of available preconditioners


    .. cpp:function:: static GenericLinearAlgebraFactory& factory()
    
        Return instance of default backend


