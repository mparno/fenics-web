
.. Documentation for the header file dolfin/mesh/SubDomain.h

.. _programmers_reference_cpp_mesh_subdomain:

SubDomain.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: SubDomain

    This class defines the interface for definition of subdomains.
    Alternatively, subdomains may be defined by a :cpp:class:`Mesh` and a
    :cpp:class:`MeshFunction` <std::size_t> over the mesh.


    .. cpp:function:: SubDomain(const double map_tol=1.0e-10)
    
        Constructor
        
        *Arguments*
            map_tol (double)
                The tolerance used when identifying mapped points using
                the function SubDomain::map.


    .. cpp:function:: bool inside(const Array<double>& x, bool on_boundary) const
    
        Return true for points inside the subdomain
        
        *Arguments*
            x (:cpp:class:`Array` <double>)
                The coordinates of the point.
            on_boundary (bool)
                True for points on the boundary.
        
        *Returns*
            bool
                True for points inside the subdomain.


    .. cpp:function:: void map(const Array<double>& x, Array<double>& y) const
    
        Map coordinate x in domain H to coordinate y in domain G (used for
        periodic boundary conditions)
        
        *Arguments*
            x (:cpp:class:`Array` <double>)
                The coordinates in domain H.
            y (:cpp:class:`Array` <double>)
                The coordinates in domain G.


    .. cpp:function:: void snap(Array<double>& x) const
    
        Snap coordinate to boundary of subdomain
        
        *Arguments*
            x (:cpp:class:`Array` <double>)
                The coordinates.


    .. cpp:function:: void mark_cells(Mesh& mesh, std::size_t sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (std::size_t) on cells for given subdomain number
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to be marked.
            sub_domain (std::size_t)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark_facets(Mesh& mesh, std::size_t sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (std::size_t) on facets for given subdomain number
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to be marked.
            sub_domain (std::size_t)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(Mesh& mesh, std::size_t dim, std::size_t sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (std::size_t) for given topological dimension
        and subdomain number
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh to be marked.
            dim (std::size_t)
                The topological dimension of entities to be marked.
            sub_domain (std::size_t)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshFunction<std::size_t>& sub_domains, std::size_t sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (std::size_t) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshFunction` <std::size_t>)
                The subdomain markers.
            sub_domain (std::size_t)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshFunction<int>& sub_domains, int sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (int) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshFunction` <int>)
                The subdomain markers.
            sub_domain (int)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshFunction<double>& sub_domains, double sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (double) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshFunction` <double>)
                The subdomain markers.
            sub_domain (double)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshFunction<bool>& sub_domains, bool sub_domain, bool check_midpoint=true) const
    
        Set subdomain markers (bool) for given subdomain
        
        *Arguments*
            sub_domains (:cpp:class:`MeshFunction` <bool>)
                The subdomain markers.
            sub_domain (bool)
                The subdomain number.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshValueCollection<std::size_t>& sub_domains, std::size_t sub_domain, const Mesh& mesh, bool check_midpoint=true) const
    
        Set subdomain markers (std::size_t) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshValueCollection` <std::size_t>)
                The subdomain markers.
            sub_domain (std::size_t)
                The subdomain number.
            mesh (:cpp:class:`Mesh`)
                The mesh.
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshValueCollection<int>& sub_domains, int sub_domain, const Mesh& mesh, bool check_midpoint=true) const
    
        Set subdomain markers (int) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshValueCollection` <int>)
                The subdomain markers
            sub_domain (int)
                The subdomain number
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshValueCollection<double>& sub_domains, double sub_domain, const Mesh& mesh, bool check_midpoint=true) const
    
        Set subdomain markers (double) for given subdomain number
        
        *Arguments*
            sub_domains (:cpp:class:`MeshValueCollection` <double>)
                The subdomain markers.
            sub_domain (double)
                The subdomain number
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: void mark(MeshValueCollection<bool>& sub_domains, bool sub_domain, const Mesh& mesh, bool check_midpoint=true) const
    
        Set subdomain markers (bool) for given subdomain
        
        *Arguments*
            sub_domains (:cpp:class:`MeshValueCollection` <bool>)
                The subdomain markers
            sub_domain (bool)
                The subdomain number
            check_midpoint (bool)
                Flag for whether midpoint of cell should be checked (default).


    .. cpp:function:: std::size_t geometric_dimension() const
    
        Return geometric dimension
        
        *Returns*
            std::size_t
                The geometric dimension.


    .. cpp:function:: void apply_markers(S& sub_domains, T sub_domain, const Mesh& mesh, bool check_midpoint) const
    
        Apply marker of type T (most likely an std::size_t) to object of class
        S (most likely MeshFunction or MeshValueCollection)


