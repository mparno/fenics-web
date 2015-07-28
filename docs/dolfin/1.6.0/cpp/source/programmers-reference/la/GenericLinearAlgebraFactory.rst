
.. Documentation for the header file dolfin/la/GenericLinearAlgebraFactory.h

.. _programmers_reference_cpp_la_genericlinearalgebrafactory:

GenericLinearAlgebraFactory.h
=============================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: GenericLinearAlgebraFactory

    .. cpp:function:: GenericLinearAlgebraFactory()
    
        Constructor


    .. cpp:function:: std::shared_ptr<GenericMatrix> create_matrix() const = 0
    
        Create empty matrix


    .. cpp:function:: std::shared_ptr<GenericVector> create_vector() const = 0
    
        Create empty vector


    .. cpp:function:: std::shared_ptr<TensorLayout> create_layout(std::size_t rank) const = 0
    
        Create empty tensor layout


    .. cpp:function:: std::shared_ptr<GenericLinearOperator> create_linear_operator() const = 0
    
        Create empty linear operator


    .. cpp:function:: std::shared_ptr<GenericLUSolver> create_lu_solver(std::string method) const = 0
    
        Create LU solver


    .. cpp:function:: std::shared_ptr<GenericLinearSolver> create_krylov_solver(std::string method, std::string preconditioner) const = 0
    
        Create Krylov solver


    .. cpp:function:: std::map<std::string, std::string> lu_solver_methods() const
    
        Return a list of available LU solver methods.  This function
        should be overloaded by subclass if non-empty.


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_methods() const
    
        Return a list of available Krylov solver methods.  This
        function should be overloaded by subclass if non-empty.


    .. cpp:function:: std::map<std::string, std::string> krylov_solver_preconditioners() const
    
        Return a list of available preconditioners.
        This function should be overloaded by subclass if non-empty.


.. cpp:class:: NotImplementedLinearOperator

    *Parent class(es)*
    
        * :cpp:class:`GenericLinearOperator`
        
