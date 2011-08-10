
.. Documentation for the header file dolfin/mesh/DomainBoundary.h

.. _programmers_reference_cpp_mesh_domainboundary:

DomainBoundary.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: DomainBoundary

    *Parent class(es)*
    
        * :cpp:class:`SubDomain`
        
    This class provides a SubDomain which picks out the boundary of
    a mesh, and provides a convenient way to specify boundary
    conditions on the entire boundary of a mesh.


    .. cpp:function:: DomainBoundary()
    
        Constructor


    .. cpp:function:: bool inside(const Array<double>& x, bool on_boundary) const
    
        Return true for points on the boundary


