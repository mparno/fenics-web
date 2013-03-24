
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
            local_mesh_data (:cpp:class:`LocalMeshData`)
                Data from which to build the mesh.


    .. cpp:function:: Mesh(const CSGGeometry& geometry, std::size_t resolution)
    
        Create mesh defined by Constructive Solid Geometry (CSG)
        
        *Arguments*
            geometry (:cpp:class:`CSGGeometry`)
                The CSG geometry
            resolution (std::size_t)
                An integer specifying the mesh resolution


    .. cpp:function:: Mesh(boost::shared_ptr<const CSGGeometry> geometry, std::size_t resolution)
    
        Create mesh defined by Constructive Solid Geometry (CSG)
        
        *Arguments*
            geometry (:cpp:class:`CSGGeometry`)
                The CSG geometry
            resolution (std::size_t)
                An integer specifying the mesh resolution


    .. cpp:function:: const Mesh& operator=(const Mesh& mesh)
    
        Assignment operator
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                Another :cpp:class:`Mesh` object.


    .. cpp:function:: std::size_t num_vertices() const
    
        Get number of vertices in mesh.
        
        *Returns*
            std::size_t
                Number of vertices.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t num_edges() const
    
        Get number of edges in mesh.
        
        *Returns*
            std::size_t
                Number of edges.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t num_faces() const
    
        Get number of faces in mesh.
        
        *Returns*
            std::size_t
                Number of faces.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t num_facets() const
    
        Get number of facets in mesh.
        
        *Returns*
            std::size_t
                Number of facets.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t num_cells() const
    
        Get number of cells in mesh.
        
        *Returns*
            std::size_t
                Number of cells.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t num_entities(std::size_t d) const
    
        Get number of entities of given topological dimension.
        
        *Arguments*
            d (std::size_t)
                Topological dimension.
        
        *Returns*
            std::size_t
                Number of entities of topological dimension d.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::vector<double>& coordinates()
    
        Get vertex coordinates.
        
        *Returns*
            std::vector<double>&
                Coordinates of all vertices.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: const std::vector<double>& coordinates() const
    
        Return coordinates of all vertices (const version).


    .. cpp:function:: const std::vector<unsigned int>& cells() const
    
        Get cell connectivity.
        
        *Returns*
            std::vector<std::size_t>
                Connectivity for all cells.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t size(std::size_t dim) const
    
        Get number of local entities of given topological dimension.
        
        *Arguments*
            dim (std::size_t)
                Topological dimension.
        
        *Returns*
            std::size_t
                Number of local entities of topological dimension d.
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t size_global(std::size_t dim) const
    
        Get global number of entities of given topological dimension.
        
        *Arguments*
            dim (std::size_t)
                Topological dimension.
        
        *Returns*
            std::size_t
                Global number of entities of topological dimension d.
        
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


    .. cpp:function:: CellType& type()
    
        Get mesh cell type.
        
        *Returns*
            :cpp:class:`CellType`
                The cell type object associated with the mesh.


    .. cpp:function:: const CellType& type() const
    
        Get mesh cell type (const version).


    .. cpp:function:: std::size_t init(std::size_t dim) const
    
        Compute entities of given topological dimension.
        
        *Arguments*
            dim (std::size_t)
                Topological dimension.
        
        *Returns*
            std::size_t
                Number of created entities.


    .. cpp:function:: void init(std::size_t d0, std::size_t d1) const
    
        Compute connectivity between given pair of dimensions.
        
        *Arguments*
            d0 (std::size_t)
                Topological dimension.
        
            d1 (std::size_t)
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


    .. cpp:function:: void rotate(double angle, std::size_t axis=2)
    
        Rotate mesh around a coordinate axis through center of mass
        of all mesh vertices
        
        *Arguments*
            angle (double)
                The number of degrees (0-360) of rotation
            axis (std::size_t)
                The coordinate axis around which to rotate the mesh


    .. cpp:function:: void rotate(double angle, std::size_t axis, const Point& p)
    
        Rotate mesh around a coordinate axis through a given point
        
        *Arguments*
            angle (double)
                The number of degrees (0-360) of rotation
            axis (std::size_t)
                The coordinate axis around which to rotate the mesh
            point (:cpp:class:`Point`)
                The point around which to rotate the mesh


    .. cpp:function:: boost::shared_ptr<MeshDisplacement> move(BoundaryMesh& boundary)
    
        Move coordinates of mesh according to new boundary coordinates.
        
        *Arguments*
            boundary (:cpp:class:`BoundaryMesh`)
                A mesh containing just the boundary cells.
        
        *Returns*
            MeshDisplacement
                Displacement encapsulated in Expression subclass MeshDisplacement


    .. cpp:function:: boost::shared_ptr<MeshDisplacement> move(Mesh& mesh)
    
        Move coordinates of mesh according to adjacent mesh with common global
        vertices.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                A :cpp:class:`Mesh` object.
        
        *Returns*
            MeshDisplacement
                Displacement encapsulated in Expression subclass MeshDisplacement


    .. cpp:function:: void move(const GenericFunction& displacement)
    
        Move coordinates of mesh according to displacement function.
        
        *Arguments*
            displacement (:cpp:class:`GenericFunction`)
                A :cpp:class:`GenericFunction` object.


    .. cpp:function:: void smooth(std::size_t num_iterations=1)
    
        Smooth internal vertices of mesh by local averaging.
        
        *Arguments*
            num_iterations (std::size_t)
                Number of iterations to perform smoothing,
                default value is 1.


    .. cpp:function:: void smooth_boundary(std::size_t num_iterations=1, bool harmonic_smoothing=true)
    
        Smooth boundary vertices of mesh by local averaging.
        
        *Arguments*
            num_iterations (std::size_t)
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


    .. cpp:function:: const std::vector<std::size_t>& color(std::string coloring_type) const
    
        Color the cells of the mesh such that no two neighboring cells
        share the same color. A colored mesh keeps a
        CellFunction<std::size_t> named "cell colors" as mesh data which
        holds the colors of the mesh.
        
        *Arguments*
            coloring_type (std::string)
                Coloring type, specifying what relation makes two
                cells neighbors, can be one of "vertex", "edge" or
                "facet".
        
        *Returns*
            MeshFunction<std::size_t>
                The colors as a mesh function over the cells of the mesh.


    .. cpp:function:: const std::vector<std::size_t>& color(std::vector<std::size_t> coloring_type) const
    
        Color the cells of the mesh such that no two neighboring cells
        share the same color. A colored mesh keeps a
        CellFunction<std::size_t> named "cell colors" as mesh data which
        holds the colors of the mesh.
        
        *Arguments*
            coloring_type (std::vector<std::size_t>)
                Coloring type given as list of topological dimensions,
                specifying what relation makes two mesh entinties neighbors.
        
        *Returns*
            MeshFunction<std::size_t>
                The colors as a mesh function over entities of the mesh.


    .. cpp:function:: void intersected_cells(const Point& point, std::set<std::size_t>& cells) const
    
        Compute all cells which are intersected by the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
            cells (std::set<std::size_t>)
                A set of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const std::vector<Point>& points, std::set<std::size_t>& cells) const
    
        Compute all cells which are intersected by any of a vector of points.
        
        *Arguments*
            points (std::vector<:cpp:class:`Point`>)
                A vector of :cpp:class:`Point` objects.
        
            cells (std::set<std::size_t>)
                A set of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const MeshEntity& entity, std::vector<std::size_t>& cells) const
    
        Compute all cells which are intersected by the given entity.
        
        *Arguments*
            entity (:cpp:class:`MeshEntity`)
                A :cpp:class:`MeshEntity` object.
        
            cells (std::vector<std::size_t>)
                A vector of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const std::vector<MeshEntity>& entities, std::set<std::size_t>& cells) const
    
        Compute all cells which are intersected by any of a vector of entities.
        
        *Arguments*
            entities (std::vector<:cpp:class:`MeshEntity`>)
                A vector of :cpp:class:`MeshEntity` objects.
        
            cells (std::set<std::size_t>)
                A vector of indices of all intersected cells.


    .. cpp:function:: void intersected_cells(const Mesh& mesh, std::set<std::size_t>& cells) const
    
        Compute all cells which are intersected by the given mesh.
        
        *Arguments*
            mesh (:cpp:class:`Mesh`)
                A :cpp:class:`Mesh` object.
        
            cells (std::set<std::size_t>)
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


    .. cpp:function:: std::size_t closest_cell(const Point& point) const
    
        Find the cell in the mesh closest to the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            std::size_t
                The index of the closest cell.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh(1, 1);
                Point point(0.0, 2.0);
                info("%d", mesh.closest_cell(point));
        
            output::
        
                1


    .. cpp:function:: std::pair<Point, std::size_t> closest_point_and_cell(const Point& point) const
    
        Find the point and corresponding cell closest to the given point.
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            std::pair<:cpp:class:`Point`, std::size_t>
                A pair consisting of the closest point and corresponding cell index.


    .. cpp:function:: double distance(const Point& point) const
    
        Computes the distance between a given point and the mesh
        
        *Arguments*
            point (:cpp:class:`Point`)
                A :cpp:class:`Point` object.
        
        *Returns*
            double
                The distance to the mesh.


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


    .. cpp:function:: double rmin() const
    
        Compute minimum cell inradius.
        
        *Returns*
            double
                The minimum of cells' inscribed sphere radii
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: double rmax() const
    
        Compute maximum cell inradius.
        
        *Returns*
            double
                The maximum of cells' inscribed sphere radii
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: double radius_ratio_min() const
    
        Compute minimum normalized radius ratio of cells.
        
        *Returns*
            double
                The minimum over cells of normalized cell
                radius ratio (which is = cell_dimension *
                * inradius / circumradius; cell_dimension
                is normalization factor).
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: double radius_ratio_max() const
    
        Compute maximum normalized radius ratio of cells.
        
        *Returns*
            double
                The maximum over cells of normalized cell
                radius ratio (which is = cell_dimension *
                * inradius / circumradius; cell_dimension
                is normalization factor).
        
        *Example*
            .. note::
        
                No example code available for this function.


    .. cpp:function:: std::size_t hash() const
    
        Compute hash of mesh, currently based on the has of the mesh
        geometry and mesh topology.
        
        *Returns*
            std::size_t
                A tree-hashed value of the coordinates over all MPI processes
        


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


    .. cpp:function:: std::vector<int>& cell_orientations()
    
        Return cell_orientations
        
        *Returns*
            std::vector<int>
                Map from cell index to orientation of cell


    .. cpp:function:: const std::vector<int>& cell_orientations() const
    
        Return cell_orientations (const version)
        
        *Returns*
            std::vector<int>
                Map from cell index to orientation of cell


    .. cpp:function:: void init_cell_orientations(const Expression& global_normal)
    
        Compute and initialize cell_orientations relative to a given
        global outward direction/normal/orientation. Only defined if
        mesh is orientable.
        
        *Arguments*
            global_normal (Expression)
                A global normal direction to the mesh


