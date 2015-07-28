
.. Documentation for the header file dolfin/mesh/MultiMesh.h

.. _programmers_reference_cpp_mesh_multimesh:

MultiMesh.h
===========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MultiMesh

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
    This class represents a collection of meshes with arbitrary
    overlaps. A multimesh may be created from a set of standard
    meshes spaces by repeatedly calling add(), followed by a call to
    build(). Note that a multimesh is not useful until build() has
    been called.


    .. cpp:function:: MultiMesh()
    
        Create empty multimesh


    .. cpp:function:: std::size_t num_parts() const
    
        Return the number of meshes (parts) of the multimesh
        
        *Returns*
            std::size_t
                The number of meshes (parts) of the multimesh.


    .. cpp:function:: std::shared_ptr<const Mesh> part(std::size_t i) const
    
        Return mesh (part) number i
        
        *Arguments*
            i (std::size_t)
                The part number
        
        *Returns*
            :cpp:class:`Mesh`
                Mesh (part) number i


    .. cpp:function:: const std::vector<unsigned int>& uncut_cells(std::size_t part) const
    
        Return the list of uncut cells for given part. The uncut cells
        are defined as all cells that don't collide with any cells in
        any other part with higher part number.
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of uncut cell indices for given part


    .. cpp:function:: const std::vector<unsigned int>& cut_cells(std::size_t part) const
    
        Return the list of cut cells for given part. The cut cells are
        defined as all cells that collide with the boundary of any
        part with higher part number.
        
        FIXME: Figure out whether this makes sense; a cell may collide
        with the boundary of part j but may still be covered
        completely by the domain of part j + 1. Possible solution is
        to for each part i check overlapping parts starting from the
        top and working back down to i + 1.
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of cut cell indices for given part


    .. cpp:function:: const std::vector<unsigned int>& covered_cells(std::size_t part) const
    
        Return the list of covered cells for given part. The covered
        cells are defined as all cells that collide with the domain of
        any part with higher part number, but not with the boundary of
        that part; in other words cells that are completely covered by
        any other part (and which therefore are inactive).
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::vector<unsigned int>
                List of covered cell indices for given part


    .. cpp:function:: const std::map<unsigned int, std::vector<std::pair<std::size_t, unsigned int> > >& collision_map_cut_cells(std::size_t part) const
    
        Return the collision map for cut cells of the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::vector<std::pair<std::size_t, unsigned int> > >
                A map from cell indices of cut cells to a list of
                cutting cells. Each cutting cell is represented as a
                pair (part_number, cutting_cell_index).


    .. cpp:function:: const std::map<unsigned int, quadrature_rule >& quadrature_rule_cut_cells(std::size_t part) const
    
        Return quadrature rules for cut cells on the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::pair<std::vector<double>, std::vector<double> > >
                A map from cell indices of cut cells to quadrature
                rules. Each quadrature rule is represented as a pair
                of a flattened array of quadrature points and a
                corresponding array of quadrature weights.


    .. cpp:function:: quadrature_rule quadrature_rule_cut_cell(std::size_t part, unsigned int cell_index) const
    
        Return quadrature rule for a given cut cell on the given part
        
        *Arguments*
            part (std::size_t)
                The part number
            cell (unsigned int)
                The cell index
        
        *Returns*
            std::pair<std::vector<double>, std::vector<double> >
                A quadrature rule represented as a pair of a flattened
                array of quadrature points and a corresponding array
                of quadrature weights. An error is raised if the given
                cell is not in the map.
        
        Developer note: this function is mainly useful from Python and
        could be replaced by a suitable typemap that would make the
        previous more general function accessible from Python.


    .. cpp:function:: const std::map<unsigned int, std::vector<quadrature_rule> >& quadrature_rule_overlap(std::size_t part) const
    
        Return quadrature rules for the overlap on the given part.
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::pair<std::vector<double>, std::vector<double> > >
                A map from cell indices of cut cells to quadrature
                rules.  A separate quadrature rule is given for each
                cutting cell and stored in the same order as in the
                collision map. Each quadrature rule is represented as
                a pair of an array of quadrature points and a
                corresponding flattened array of quadrature weights.


    .. cpp:function:: const std::map<unsigned int, std::vector<quadrature_rule> >& quadrature_rule_interface(std::size_t part) const
    
        Return quadrature rules for the interface on the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::pair<std::vector<double>, std::vector<double> > >
                A map from cell indices of cut cells to quadrature
                rules on an interface part cutting through the cell.
                A separate quadrature rule is given for each cutting
                cell and stored in the same order as in the collision
                map. Each quadrature rule is represented as a pair of
                an array of quadrature points and a corresponding
                flattened array of quadrature weights.


    .. cpp:function:: const std::map<unsigned int, std::vector<std::vector<double> > >& facet_normals(std::size_t part) const
    
        Return facet normals for the interface on the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::map<unsigned int, std::vector<std::vector<double> > >
                A map from cell indices of cut cells to facet normals
                on an interface part cutting through the cell. A
                separate list of facet normals, one for each
                quadrature point, is given for each cutting cell and
                stored in the same order as in the collision map. The
                facet normals for each set of quadrature points is
                stored as a contiguous flattened array, the length of
                which should be equal to the number of quadrature
                points multiplied by the geometric dimension. Puh!


    .. cpp:function:: std::shared_ptr<const BoundingBoxTree> bounding_box_tree(std::size_t part) const
    
        Return the bounding box tree for the mesh of the given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::shared_ptr<const BoundingBoxTree>
                The bounding box tree


    .. cpp:function:: std::shared_ptr<const BoundingBoxTree> bounding_box_tree_boundary(std::size_t part) const
    
        Return the bounding box tree for the boundary mesh of the
        given part
        
        *Arguments*
            part (std::size_t)
                The part number
        
        *Returns*
            std::shared_ptr<const BoundingBoxTree>
                The bounding box tree


    .. cpp:function:: void add(std::shared_ptr<const Mesh> mesh)
    
        Add mesh (shared pointer version)
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh


    .. cpp:function:: void add(const Mesh& mesh)
    
        Add mesh (reference version)
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh


    .. cpp:function:: void build()
    
        Build multimesh


    .. cpp:function:: void clear()
    
        Clear multimesh


    .. cpp:function:: static Parameters default_parameters()
    
        Default parameter values


