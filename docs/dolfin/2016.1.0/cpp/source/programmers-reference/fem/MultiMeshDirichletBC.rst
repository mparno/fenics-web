
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


    .. cpp:function:: MultiMeshDirichletBC(std::shared_ptr<const MultiMeshFunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::shared_ptr<const SubDomain> sub_domain, std::string method="topological", bool check_midpoint=true, bool exclude_overlapped_boundaries=true)
    
        Create boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`MultiMeshFunctionSpace`)
                The function space
            g (:cpp:class:`GenericFunction`)
                The value
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain
            method (std::string)
                Option passed to DirichletBC.
            check_midpoint (bool)
                Option passed to DirichletBC.
            exclude_overlapped_boundaries (bool)
                If true, then the variable on_boundary will
                be set to false for facets that are overlapped
                by another mesh (irrespective of the layering order
                of the meshes).


    .. cpp:function:: MultiMeshDirichletBC(std::shared_ptr<const MultiMeshFunctionSpace> V, std::shared_ptr<const GenericFunction> g, std::shared_ptr<const MeshFunction<std::size_t>> sub_domains, std::size_t sub_domain, std::size_t part, std::string method="topological")
    
        Create boundary condition for subdomain specified by index
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            g (:cpp:class:`GenericFunction`)
                The value.
            sub_domains (:cpp:class:`MeshFunction` <std::size_t>)
                Subdomain markers
            sub_domain (std::size_t)
                The subdomain index (number)
            part (std::size_t)
                The part on which to set boundary conditions
            method (std::string)
                Optional argument: A string specifying the
                method to identify dofs.


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
        
