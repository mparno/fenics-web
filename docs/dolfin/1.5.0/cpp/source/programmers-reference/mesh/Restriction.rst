
.. Documentation for the header file dolfin/mesh/Restriction.h

.. _programmers_reference_cpp_mesh_restriction:

Restriction.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Restriction

    This class represents a restriction of a mesh to a subdomain,
    which can be defined as a subset of all the cells, the facets,
    or possibly lower dimensional entities of the mesh.


    .. cpp:function:: Restriction(const Mesh& mesh, const SubDomain& sub_domain)
    
        Create cell-based restriction from subdomain
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh
            sub_domain (:cpp:class:`SubDomain`)
                Sub domain defining the restriction


    .. cpp:function:: Restriction(const Mesh& mesh, const SubDomain& sub_domain, std::size_t dim)
    
        Create restriction from subdomain to entities of arbitrary dimension
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh
            sub_domain (:cpp:class:`SubDomain`)
                Sub domain defining the restriction
            dim (std::size_t)
                Dimension of restriction


    .. cpp:function:: Restriction(const MeshFunction<std::size_t>& domain_markers, std::size_t domain_number)
    
        Create restriction from domain markers
        
        *Arguments*
            domain_markers (:cpp:class:`MeshFunction` <std::size_t>)
                Domain markers for the cells of the mesh.
            domain_number (std::size_t)
                Identifier for domain.


    .. cpp:function:: Restriction(std::shared_ptr<const MeshFunction<std::size_t> > domain_markers, std::size_t domain_number)
    
        Create restriction from domain markers (shared pointer version)
        
        *Arguments*
            domain_markers (:cpp:class:`MeshFunction` <std::size_t>)
                Domain markers for the cells of the mesh.
            domain_number (std::size_t)
                Identifier for domain.


    .. cpp:function:: const Mesh& mesh() const
    
        Return the full unrestricted mesh


    .. cpp:function:: std::size_t dim() const
    
        Return topological dimension of restriction


    .. cpp:function:: bool contains(const MeshEntity& entity) const
    
        Check whether restriction contains entity


    .. cpp:function:: bool contains(std::size_t d, std::size_t i) const
    
        Check whether restriction contains entity (d, i)


