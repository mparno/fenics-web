
.. Documentation for the header file dolfin/fem/EqualityBC.h

.. _programmers_reference_cpp_fem_equalitybc:

EqualityBC.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EqualityBC

    *Parent class(es)*
    
        * :cpp:class:`BoundaryCondition`
        
    This class specifies the interface for setting equality boundary
    conditions for partial differential equations,
    
       u(x) = u(y),    for all x and y on G,
    
    where G is subdomain of the mesh.
    
    The sub domain G may be specified in two different ways. Both of
    them produce a set of unknowns (dofs) with should be equal.
    
    The simplest approach is to specify a SubDomain object, using
    the inside() function to specify on which facets the boundary
    condition should be applied.
    
    Alternatively, the boundary may be specified by the boundary
    indicators included in the mesh.
    
    Current implementation assume that the problem is scalar,
    so in case of mixed systems (vector-valued and mixed elements)
    all compoments will be set equal.


