
.. Documentation for the header file dolfin/mesh/Facet.h

.. _programmers_reference_cpp_mesh_facet:

Facet.h
=======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Facet

    *Parent class(es)*
    
        * :cpp:class:`MeshEntity`
        
    A Facet is a MeshEntity of topological codimension 1.


    .. cpp:function:: Facet(const Mesh& mesh, std::size_t index)
    
        Constructor


    .. cpp:function:: double normal(std::size_t i) const
    
        Compute component i of the normal to the facet


    .. cpp:function:: Point normal() const
    
        Compute normal to the facet


    .. cpp:function:: bool exterior() const
    
        Return true if facet is an exterior facet (relative to global mesh,
        so this function will return false for facets on partition
        boundaries). Facet connectivity must be initialized before
        calling this function.


    .. cpp:function:: std::pair<const Cell, const Cell> adjacent_cells(const std::vector<std::size_t>* facet_orientation) const
    
        Return adjacent cells. An optional argument that lists for
        each facet the index of the first cell may be given to specify
        the ordering of the two cells. If not specified, the ordering
        will depend on the (arbitrary) ordering of the mesh
        connectivity.


.. cpp:class:: FacetFunction

    *Parent class(es)*
    
        * :cpp:class:`MeshFunction<T>`
        
    A FacetFunction is a MeshFunction of topological codimension 1.


