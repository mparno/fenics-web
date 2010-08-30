.. Documentation for the header file dolfin/mesh/Facet.h

.. _programmers_reference_cpp_mesh_facet:

Facet.h
=======

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: Facet

    *Parent class*
    
        * :cpp:class:`MeshEntity`
        
    A Facet is a MeshEntity of topological codimension 1.

    .. cpp:function:: Facet(const Mesh& mesh, uint index)
    
        Constructor

    .. cpp:function:: bool interior() const
    
        Determine whether or not facet is an interior facet. This is 'relative'
        to the given partition of the mesh if the mesh is distributed

    .. cpp:function:: std::pair<const Cell, const Cell> adjacent_cells(MeshFunction<uint>* facet_orientation=0) const
    
        Return adjacent cells. An optional argument that lists for
        each facet the index of the first cell may be given to specify
        the ordering of the two cells. If not specified, the ordering
        will depend on the (arbitrary) ordering of the mesh
        connectivity.

.. cpp:class:: FacetIterator

    *Parent class*
    
        * :cpp:class:`MeshEntityIterator`
        
    A FacetIterator is a MeshEntityIterator of topological codimension 1.

.. cpp:class:: T>

    *Parent class*
    
        * :cpp:class:`MeshFunction<T>`
        
    A FacetFunction is a MeshFunction of topological codimension 1.

