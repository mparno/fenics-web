.. Documentation for the header file dolfin/mesh/MeshEntity.h

.. _programmers_reference_cpp_mesh_Mesh:

MeshEntity.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: MeshEntity

        A MeshEntity represents a mesh entity associated with
        a specific topological dimension of some mesh.

    .. cpp:function:: MeshEntity() :_mesh(0), _dim(0), _index(0)
    
        Default Constructor

    .. cpp:function:: MeshEntity(const Mesh& mesh, uint dim, uint index)
    
        Constructor

    .. cpp:function:: Point midpoint() const
    
        Compute midpoint of cell

    .. cpp:function:: bool incident(const MeshEntity& entity) const
    
        Check if given entity is indicent

    .. cpp:function:: bool intersects(const MeshEntity& entity) const
    
        Check if given entity intersects (using inexact but fast numerics)

    .. cpp:function:: bool intersects(const Point& point) const
    
        Check if given point intersects (using inexact but fast numerics)

    .. cpp:function:: bool intersects_exactly(const MeshEntity& entity) const
    
        Check if given entity intersects (using exact numerics)

    .. cpp:function:: bool intersects_exactly(const Point& point) const
    
        Check if given point intersects (using exact numerics)

    .. cpp:function:: bool operator==(const MeshEntity& another) const
    
        Comparision Operator

    .. cpp:function:: const Mesh& mesh() const
    
        Return mesh associated with mesh entity

    .. cpp:function:: const uint* entities(uint dim) const
    
        Return array of indices for incident mesh entitites of given topological dimension

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

    .. cpp:function:: template <typename K> CGAL::Bbox_3 bbox() const
    
        Returns a 3D bounding box of the mesh entity. For lower dimension it may be a degenerated box.

    .. cpp:function:: uint dim() const
    
        Return topological dimension

    .. cpp:function:: uint index() const
    
        Return index of mesh entity

    .. cpp:function:: uint index(const MeshEntity& entity) const
    
        Compute local index of given incident entity (error if not found)

    .. cpp:function:: uint num_entities(uint dim) const
    
        Return number of incident mesh entities of given topological dimension

    .. cpp:function:: virtual ~MeshEntity()
    
        Destructor

