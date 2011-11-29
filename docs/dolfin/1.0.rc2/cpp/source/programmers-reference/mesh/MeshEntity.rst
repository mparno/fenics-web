
.. Documentation for the header file dolfin/mesh/MeshEntity.h

.. _programmers_reference_cpp_mesh_meshentity:

MeshEntity.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshEntity

    A MeshEntity represents a mesh entity associated with
    a specific topological dimension of some :cpp:class:`Mesh`.


    .. cpp:function:: MeshEntity()
    
        Default Constructor


    .. cpp:function:: MeshEntity(const Mesh& mesh, uint dim, uint index)
    
        Constructor
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            dim (uint)
                The topological dimension.
            index (uint)
                The index.


    .. cpp:function:: void init(const Mesh& mesh, uint dim, uint index)
    
        Initialize mesh entity with given data
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                The mesh.
            dim (uint)
                The topological dimension.
            index (uint)
                The index.


    .. cpp:function:: bool operator==(const MeshEntity& another) const
    
        Comparision Operator
        
        *Arguments*
            another (:cpp:class:`MeshEntity`)
                Another mesh entity
        
        *Returns*
            bool
                True if the two mesh entities are equal.


    .. cpp:function:: bool operator!=(const MeshEntity& another) const
    
        Comparision Operator
        
        *Arguments*
            another (MeshEntity)
                Another mesh entity.
        
        *Returns*
            bool
                True if the two mesh entities are NOT equal.


    .. cpp:function:: const Mesh& mesh() const
    
        Return mesh associated with mesh entity
        
        *Returns*
            :cpp:class:`Mesh`
                The mesh.


    .. cpp:function:: uint dim() const
    
        Return topological dimension
        
        *Returns*
            uint
                The dimension.


    .. cpp:function:: uint index() const
    
        Return index of mesh entity
        
        *Returns*
            uint
                The index.


    .. cpp:function:: uint num_entities(uint dim) const
    
        Return number of incident mesh entities of given topological dimension
        
        *Arguments*
            dim (uint)
                The topological dimension.
        
        *Returns*
            uint
                The number of incident MeshEntity objects of given dimension.


    .. cpp:function:: const uint* entities(uint dim) const
    
        Return array of indices for incident mesh entitites of given
        topological dimension
        
        *Arguments*
            dim (uint)
                The topological dimension.
        
        *Returns*
            uint
                The index for incident mesh entities of given dimension.


    .. cpp:function:: uint mesh_id() const
    
        Return unique mesh ID
        
        *Returns*
            uint
                The unique mesh ID.


    .. cpp:function:: bool incident(const MeshEntity& entity) const
    
        Check if given entity is incident
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The entity.
        
        *Returns*
            bool
                True if the given entity is incident


    .. cpp:function:: bool intersects(const Point& point) const
    
        Check if given point intersects (using inexact but fast
        numerics)
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
        
        *Returns*
            bool
                True if the given point intersects.


    .. cpp:function:: bool intersects(const MeshEntity& entity) const
    
        Check if given entity intersects (using inexact but fast
        numerics)
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The mesh entity.
        
        *Returns*
            bool
                True if the given entity intersects.


    .. cpp:function:: bool intersects_exactly(const Point& point) const
    
        Check if given point intersects (using exact numerics)
        
        *Arguments*
            point (:cpp:class:`Point`)
                The point.
        
        *Returns*
            bool
                True if the given point intersects.


    .. cpp:function:: bool intersects_exactly(const MeshEntity& entity) const
    
        Check if given entity intersects (using exact numerics)
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The mesh entity.
        
        *Returns*
            bool
                True if the given entity intersects.


    .. cpp:function:: uint index(const MeshEntity& entity) const
    
        Compute local index of given incident entity (error if not
        found)
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                The mesh entity.
        
        *Returns*
            uint
                The local index of given entity.


    .. cpp:function:: Point midpoint() const
    
        Compute midpoint of cell
        
        *Returns*
            :cpp:class:`Point`
                The midpoint of the cell.


    .. cpp:function:: CGAL::Bbox_3 bbox() const
    
        Returns a 3D bounding box of the mesh entity. For lower
        dimension it may be a degenerated box.


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the function space.


