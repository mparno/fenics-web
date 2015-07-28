
.. Documentation for the header file dolfin/fem/MultiMeshDirichletBC.h

.. _programmers_reference_cpp_fem_multimeshdirichletbc:

MultiMeshDirichletBC.h
======================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMeshDirichletBC

    This class is used to set Dirichlet boundary conditions for
    multimesh function spaces.


    .. cpp:function:: MultiMeshDirichletBC(const MultiMeshFunctionSpace& V, const GenericFunction& g, const SubDomain& sub_domain, std::string method="topological", bool check_midpoint=true)
    
        Create boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain.
            method (std::string)
                Optional argument: A string specifying
                the method to identify dofs.


    .. cpp:function:: MultiMeshDirichletBC(std::shared_ptr<const MultiMeshFunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::shared_ptr<const SubDomain> sub_domain, std::string method="topological", bool check_midpoint=true)
    
        Create boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The function space
            g (:cpp:class:`GenericFunction`)
                The value
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain
            method (std::string)
                Optional argument: A string specifying
                the method to identify dofs


    .. cpp:function:: void apply(GenericMatrix& A) const
    
        Apply boundary condition to a matrix
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary condition to.


    .. cpp:function:: void apply(GenericVector& b) const
    
        Apply boundary condition to a vector
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary condition to.


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b) const
    
        Apply boundary condition to a linear system
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary condition to.
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary condition to.


    .. cpp:function:: void apply(GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to vectors for a nonlinear problem
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary conditions to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a linear system for a nonlinear problem
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply boundary conditions to.
            b (:cpp:class:`GenericVector`)
                The vector to apply boundary conditions to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


.. cpp:class:: MultiMeshSubDomain

    *Parent class(es)*
    
        * :cpp:class:`SubDomain`
        
