
.. Documentation for the header file dolfin/mesh/FacetCell.h

.. _programmers_reference_cpp_mesh_facetcell:

FacetCell.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: FacetCell

    *Parent class(es)*
    
        * :cpp:class:`Cell`
        
    This class represents a cell in a mesh incident to a facet on
    the boundary. It is useful in cases where one needs to iterate
    over a boundary mesh and access the corresponding cells in the
    original mesh.


    .. cpp:function:: FacetCell(const BoundaryMesh& mesh, const Cell& facet)
    
        Create cell on mesh corresponding to given facet (cell) on boundary


    .. cpp:function:: std::size_t facet_index() const
    
        Return local index of facet with respect to the cell


