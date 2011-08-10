
.. Documentation for the header file dolfin/fem/PeriodicBC.h

.. _programmers_reference_cpp_fem_periodicbc:

PeriodicBC.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: PeriodicBC

    *Parent class(es)*
    
        * :cpp:class:`BoundaryCondition`
        
    This class specifies the interface for setting periodic boundary
    conditions for partial differential equations,
    
    .. math::
    
        u(x) &= u(F^{-1}(x)) \hbox { on } G,
    
        u(x) &= u(F(x))      \hbox{ on } H,
    
    where F : H --> G is a map from a subdomain H to a subdomain G.
    
    A periodic boundary condition must be defined by the domain G
    and the map F pulling coordinates back from H to G. The domain
    and the map are both defined by a subclass of :cpp:class:`SubDomain` which
    must overload both the inside() function, which specifies the
    points of G, and the map() function, which specifies the map
    from the points of H to the points of G.
    
    The implementation is based on matching degrees of freedom on G
    with degrees of freedom on H and only works when the mapping F
    is bijective between the sets of coordinates associated with the
    two domains. In other words, the nodes (degrees of freedom) must
    be aligned on G and H.
    
    The matching of degrees of freedom is done at the construction
    of the periodic boundary condition and is reused on subsequent
    applications to a linear system. The matching may be recomputed
    by calling the ``rebuild()`` function.


    .. cpp:function:: PeriodicBC(const FunctionSpace& V, const SubDomain& sub_domain)
    
        Create periodic boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            sub_domain (:cpp:class:`SubDomain`)
                The sub domain.


    .. cpp:function:: PeriodicBC(boost::shared_ptr<const FunctionSpace> V, boost::shared_ptr<const SubDomain> sub_domain)
    
        Create periodic boundary condition for subdomain
        
        *Arguments*
            V (:cpp:class:`FunctionSpace`)
                The function space.
            sub_domain (:cpp:class:`SubDomain`)
                The subdomain.


    .. cpp:function:: void apply(GenericMatrix& A) const
    
        Apply boundary condition to a matrix
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply bc to.


    .. cpp:function:: void apply(GenericVector& b) const
    
        Apply boundary condition to a vector
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply bc to.


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b) const
    
        Apply boundary condition to a linear system
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix.
            b (:cpp:class:`GenericVector`)
                The vector.


    .. cpp:function:: void apply(GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a vector for a nonlinear problem
        
        *Arguments*
            b (:cpp:class:`GenericVector`)
                The vector to apply bc to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


    .. cpp:function:: void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const
    
        Apply boundary condition to a linear system for a nonlinear
        problem
        
        *Arguments*
            A (:cpp:class:`GenericMatrix`)
                The matrix to apply bc to.
            b (:cpp:class:`GenericVector`)
                The vector to apply bc to.
            x (:cpp:class:`GenericVector`)
                Another vector (nonlinear problem).


    .. cpp:function:: void rebuild()
    
        Rebuild mapping between dofs


