
.. Documentation for the header file dolfin/mesh/MeshDomains.h

.. _programmers_reference_cpp_mesh_meshdomains:

MeshDomains.h
=============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshDomains

    The class :cpp:class:`MeshDomains` stores the division of a :cpp:class:`Mesh` into
    subdomains. For each topological dimension 0 <= d <= D, where D
    is the topological dimension of the :cpp:class:`Mesh`, a set of integer
    markers are stored for a subset of the entities of dimension d,
    indicating for each entity in the subset the number of the
    subdomain. It should be noted that the subset does not need to
    contain all entities of any given dimension; entities not
    contained in the subset are "unmarked".


    .. cpp:function:: MeshDomains()
    
        Create empty mesh domains


    .. cpp:function:: uint dim() const
    
        Return maximal topological dimension of stored markers


    .. cpp:function:: uint num_marked(uint dim) const
    
        Return number of marked entities of given dimension


    .. cpp:function:: bool is_empty() const
    
        Check whether domain data is empty


    .. cpp:function:: MeshValueCollection<unsigned int>& markers(uint dim)
    
        Get subdomain markers for given dimension


    .. cpp:function:: const MeshValueCollection<unsigned int>& markers(uint dim) const
    
        Get subdomain markers for given dimension (const version)


    .. cpp:function:: boost::shared_ptr<MeshValueCollection<unsigned int> > markers_shared_ptr(uint dim)
    
        Get subdomain markers for given dimension (shared pointer version)


    .. cpp:function:: boost::shared_ptr<const MeshValueCollection<unsigned int> > markers_shared_ptr(uint dim) const
    
        Get subdomain markers for given dimension (const shared pointer version)


    .. cpp:function:: boost::shared_ptr<const MeshFunction<unsigned int> > cell_domains(const Mesh& mesh) const
    
        Get cell domains. This function computes the mesh function
        corresponding to markers of dimension D. The mesh function is
        cached for later access and will be computed on the first call
        to this function.


    .. cpp:function:: boost::shared_ptr<const MeshFunction<unsigned int> > facet_domains(const Mesh& mesh) const
    
        Get facet domains. This function computes the mesh function
        corresponding to markers of dimension D-1. The mesh function
        is cached for later access and will be computed on the first
        call to this function.


    .. cpp:function:: void init(uint dim)
    
        Initialize mesh domains for given topological dimension


    .. cpp:function:: void clear()
    
        Clear all data


