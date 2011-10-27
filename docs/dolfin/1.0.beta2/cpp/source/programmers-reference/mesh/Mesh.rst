
.. Documentation for the header file dolfin/mesh/Mesh.h

.. _programmers_reference_cpp_mesh_mesh:

Mesh.h
======

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: Mesh

    *Parent class(es)*
    
        * :cpp:class:`Variable`
        
        * :cpp:class:`Hierarchical<Mesh>`
        
    A :cpp:class:`Mesh` consists of a set of connected and numbered mesh entities.
    
    Both the representation and the interface are
    dimension-independent, but a concrete interface is also provided
    for standard named mesh entities:
    
    .. tabularcolumns:: |c|c|c|
    
    +--------+-----------+-------------+
    | Entity | Dimension | Codimension |
    +========+===========+=============+
    | Vertex |  0        |             |
    +--------+-----------+-------------+
    | Edge   |  1        |             |
    +--------+-----------+-------------+
    | Face   |  2        |             |
    +--------+-----------+-------------+
    | Facet  |           |      1      |
    +--------+-----------+-------------+
    | Cell   |           |      0      |
    +--------+-----------+-------------+
    
    When working with mesh iterators, all entities and connectivity
    are precomputed automatically the first time an iterator is
    created over any given topological dimension or connectivity.
    
    Note that for efficiency, only entities of dimension zero
    (vertices) and entities of the maximal dimension (cells) exist
    when creating a :cpp:class:`Mesh`. Other entities must be explicitly created
    by calling init(). For example, all edges in a mesh may be
    created by a call to mesh.init(1). Similarly, connectivities
    such as all edges connected to a given vertex must also be
    explicitly created (in this case by a call to mesh.init(0, 1)).


    .. cpp:function:: Mesh()
    
        Create empty mesh


    .. cpp:function:: Mesh(const Mesh& mesh)
    
        Copy constructor.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Object to be copied.


    .. cpp:function:: explicit Mesh(std::string filename)
    
        Create mesh from data file.
        
        *Arguments*
            filename (std::string)
                Name of file to load.


    .. cpp:function:: explicit Mesh(LocalMeshData& local_mesh_data)
    
        Create a distributed mesh from local (per process) data.
        
        *Arguments*
            local_mesh_data (LocalMeshData)
                Data from which to build the mesh.


    .. cpp:function:: const Mesh& operator=(const Mesh& mesh)
    
        Assignment operator
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Another :cpp:class:`Mesh` object.


    .. cpp:function:: uint num_vertices() const
    
        Get number of vertices in mesh.
        
        *Returns*
            uint
                Number of vertices.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint num_edges() const
    
        Get number of edges in mesh.
        
        *Returns*
            uint
                Number of edges.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint num_faces() const
    
        Get number of faces in mesh.
        
        *Returns*
            uint
                Number of faces.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint num_facets() const
    
        Get number of facets in mesh.
        
        *Returns*
            uint
                Number of facets.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint num_cells() const
    
        Get number of cells in mesh.
        
        *Returns*
            uint
                Number of cells.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint num_entities(uint d) const
    
        Get number of entities of given topological dimension.
        
        *Arguments*
            d (uint)
                Topological dimension.
        
        *Returns*
            uint
                Number of entities of topological dimension d.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: double* coordinates()
    
        Get vertex coordinates.
        
        *Returns*
            double*
                Coordinates of all vertices.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: const double* coordinates() const
    
        Return coordinates of all vertices (const version).


    .. cpp:function:: const uint* cells() const
    
        Get cell connectivity.
        
        *Returns*
            uint*
                Connectivity for all cells.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: uint size(uint dim) const
    
        Get number of entities of given topological dimension.
        
        *Arguments*
            dim (uint)
                Topological dimension.
        
        *Returns*
            uint
                Number of entities of topological dimension d.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: MeshTopology& topology()
    
        Get mesh topology.
        
        *Returns*
            :cpp:class:`MeshTopology`
                The topology object associated with the mesh.


    .. cpp:function:: const MeshTopology& topology() const
    
        Get mesh topology (const version).


    .. cpp:function:: MeshGeometry& geometry()
    
        Get mesh geometry.
        
        *Returns*
            :cpp:class:`MeshGeometry`
                The geometry object associated with the mesh.


    .. cpp:function:: const MeshGeometry& geometry() const
    
        Get mesh geometry (const version).


    .. cpp:function:: MeshDomains& domains()
    
        Get mesh (sub)domains.
        
        *Returns*
            :cpp:class:`MeshDomains`
                The (sub)domains associated with the mesh.


    .. cpp:function:: const MeshDomains& domains() const
    
        Get mesh (sub)domains.


    .. cpp:function:: IntersectionOperator& intersection_operator()
    
        Get intersection operator.
        
        *Returns*
            :cpp:class:`IntersectionOperator`
                The intersection operator object associated with the mesh.


    .. cpp:function:: const IntersectionOperator& intersection_operator() const
    
        Return intersection operator (const version);


    .. cpp:function:: MeshData& data()
    
        Get mesh data.
        
        *Returns*
            :cpp:class:`MeshData`
                The mesh data object associated with the mesh.


    .. cpp:function:: const MeshData& data() const
    
        Get mesh data (const version).


    .. cpp:function:: ParallelData& parallel_data()
    
        Get parallel mesh data.
        
        *Returns*
            :cpp:class:`ParallelData`
                The parallel data object associated with the mesh.


    .. cpp:function:: const ParallelData& parallel_data() const
    
        Get parallel mesh data (const version).


    .. cpp:function:: CellType& type()
    
        Get mesh cell type.
        
        *Returns*
            :cpp:class:`CellType`
                The cell type object associated with the mesh.


    .. cpp:function:: const CellType& type() const
    
        Get mesh cell type (const version).


    .. cpp:function:: uint init(uint dim) const
    
        Compute entities of given topological dimension.
        
        *Arguments*
            dim (uint)
                Topological dimension.
        
        *Returns*
            uint
                Number of created entities.


    .. cpp:function:: void init(uint d0, uint d1) const
    
        Compute connectivity between given pair of dimensions.
        
        *Arguments*
            d0 (uint)
                Topological dimension.
        
            d1 (uint)
                Topological dimension.


    .. cpp:function:: void init() const
    
        Compute all entities and connectivity.


    .. cpp:function:: void clear()
    
        Clear all mesh data.


    .. cpp:function:: void clean()
    
        Clean out all auxiliary topology data. This clears all
        topological data, except the connectivity between cells and
        vertices.


    .. cpp:function:: void order()
    
        Order all mesh entities.
        
        .. seealso::
        
            UFC documentation (put link here!)


    .. cpp:function:: bool ordered() const
    
        Check if mesh is ordered according to the UFC numbering convention.
        
        *Returns*
            bool
                The return values is true iff the mesh is ordered.


    .. cpp:function:: void move(BoundaryMesh& boundary)
    
        Move coordinates of mesh according to new boundary coordinates.
        
        *Arguments*
            boundary (:cpp:class:`BoundaryMesh`)
                A mesh containing just the boundary cells.


    .. cpp:function:: void move(Mesh& mesh)
    
        Move coordinates of mesh according to adjacent mesh with common global
        vertices.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                A :cpp:class:`Mesh` object.


    .. cpp:function:: void move(const Function& displacement)
    
        Move coordinates of mesh according to displacement function.
        
        *Arguments*
            displacement (:cpp:class:`Function`)
                A :cpp:class:`Function` object.


    .. cpp:function:: void smooth(uint num_iterations=1)
    
        Smooth internal vertices of mesh by local averaging.
        
        *Arguments*
            num_iterations (uint)
                Number of iterations to perform smoothing,
                default value is 1.


    .. cpp:function:: void smooth_boundary(uint num_iterations=1, bool harmonic_smoothing=true)
    
        Smooth boundary vertices of mesh by local averaging.
        
        *Arguments*
            num_iterations (uint)
                Number of iterations to perform smoothing,
                default value is 1.
        
            harmonic_smoothing (bool)
                Flag to turn on harmonics smoothing, default
                value is true.


    .. cpp:function:: void snap_boundary(const SubDomain& sub_domain, bool harmonic_smoothing=true)
    
        Snap boundary vertices of mesh to match given sub domain.
        
        *Arguments*
            sub_domain (:cpp:class:`SubDomain`)
                A :cpp:class:`SubDomain` object.
        
            harmonic_smoothing (bool)
                Flag to turn on harmonics smoothing, default
                value is true.


    .. cpp:function:: const MeshFunction<unsigned int>& color(std::string coloring_type) const
    
        Color the cells of the mesh such that no two neighboring cells
        share the same color. A colored mesh keeps a
        CellFunction<unsigned int> named "cell colors" as mesh data which
        holds the colors of the mesh.
        
        *Arguments*
            coloring_type (std::string)
                Coloring type, specifying what relation makes two
                cells neighbors, can be one of "vertex", "edge" or
                "facet".
        
        *Returns*
            MeshFunction<unsigned int>
                The colors as a mesh function over the cells of the mesh.


    .. cpp:function:: const MeshFunction<unsigned int>& color(std::vector<unsigned int> coloring_type) const
    
        Color the cells of the mesh such that no two neighboring cells
        share the same color. A colored mesh keeps a
        CellFunction<unsigned int> named "cell colors" as mesh data which
        holds the colors of the mesh.
        
        *Arguments*
            coloring_type (std::vector<unsigned int>)
                Coloring type given as list of topological dimensions,
                specifying what relation makes two mesh entinties neighbors.
        
        *Returns*
            MeshFunction<unsigned int>
                The colors as a mesh function over entities of the mesh.


    .. cpp:function:: void intersected_cells(const Point& point, std::set<uint>& cells) const
    
        Compute all cells which are intersected by the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
            cells (std::set<uint>)
                A set of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const std::vector<Point>& points, std::set<uint>& cells) const
    
        Compute all cells which are intersected by any of a vector of points.
        
        *Arguments*
            points (std::vector<:cpp:class:`Point`>)
                A vector of :cpp:class:`Point` objects.
        
            cells (std::set<uint>)
                A set of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const MeshEntity& entity, std::vector<uint>& cells) const
    
        Compute all cells which are intersected by the given entity.
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                A :cpp:class:`MeshEntity` object.
        
            cells (std::vector<uint>)
                A vector of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const std::vector<MeshEntity>& entities, std::set<uint>& cells) const
    
        Compute all cells which are intersected by any of a vector of entities.
        
        *Arguments*
            entities (std::vector<:cpp:class:`MeshEntity`>)
                A vector of :cpp:class:`MeshEntity` objects.
        
            cells (std::set<uint>)
                A vector of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const Mesh& mesh, std::set<uint>& cells) const
    
        Compute all cells which are intersected by the given mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                A :cpp:class:`Mesh` object.
        
            cells (std::set<uint>)
                A set of indices of all intersected cells.


    .. cpp:function:: int intersected_cell(const Point& point) const
    
        Find the cell (if any) containing the given point. If the point
        is contained in several cells, the first cell is returned.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            int
                The index of the cell containing the point. If no cell
                is found, the return value is -1.


    .. cpp:function:: Point closest_point(const Point& point) const
    
        Find the point in the mesh closest to the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            :cpp:class:`Point`
                The closest point.


    .. cpp:function:: dolfin::uint closest_cell(const Point& point) const
    
        Find the cell in the mesh closest to the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            uint
                The index of the closest cell.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Point point(0.0, 2.0);
                info("%d", mesh.closest_cell(point));
        
            output::
        
                1


    .. cpp:function:: std::pair<Point, dolfin::uint> closest_point_and_cell(const Point& point) const
    
        Find the point and corresponding cell closest to the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            std::pair<:cpp:class:`Point`, uint>
                A pair consisting of the closest point and corresponding cell index.


    .. cpp:function:: double hmin() const
    
        Compute minimum cell diameter.
        
        *Returns*
            double
                The minimum cell diameter, the diameter is computed as
                two times the circumradius
                (http://mathworld.wolfram.com).
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: double hmax() const
    
        Compute maximum cell diameter.
        
        *Returns*
            double
                The maximum cell diameter, the diameter is computed as
                two times the circumradius
                (http://mathworld.wolfram.com).
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::string str(bool verbose) const
    
        Informal string representation.
        
        *Arguments*
            verbose (bool)
                Flag to turn on additional output.
        
        *Returns*
            std::string
                An informal representation of the mesh.
        
        *Example*
            .. note::
        
                No example code available for this function.


