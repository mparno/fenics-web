.. Documentation for the header file dolfin/mesh/SubDomain.h

.. _programmers_reference_cpp_mesh_subdomain:

SubDomain.h
===========

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: SubDomain

        This class defines the interface for definition of sub domains.
        Alternatively, sub domains may be defined by a Mesh and a
        MeshFunction<uint> over the mesh.

    .. cpp:function:: SubDomain()
    
        Constructor

    .. cpp:function:: uint geometric_dimension() const
    
        Return geometric dimension

    .. cpp:function:: virtual bool inside(const Array<double>& x, bool on_boundary) const
    
        Return true for points inside the subdomain

    .. cpp:function:: virtual void map(const Array<double>& x, Array<double>&) const
    
        Map coordinate x in domain H to coordinate y in domain G (used for
        periodic boundary conditions)

    .. cpp:function:: virtual void snap(Array<double>& x) const
    
        Snap coordinate to boundary of sub domain

    .. cpp:function:: virtual ~SubDomain()
    
        Destructor

    .. cpp:function:: void mark(MeshFunction<uint>& sub_domains, uint sub_domain) const
    
        Set sub domain markers for given subdomain

