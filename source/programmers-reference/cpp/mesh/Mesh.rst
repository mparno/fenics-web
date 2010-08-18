.. Documentation for the header file dolfin/mesh/Mesh.h

.. _programmers_reference_cpp_mesh_Mesh:

Mesh.h
======

.. cpp:class:: Mesh

    *Parent class*

        * :cpp:class:`Variable`

    A Mesh consists of a set of connected and numbered mesh entities.

    Both the representation and the interface are dimension-independent,
    but a concrete interface is also provided for standard named mesh
    entities:

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
    | Cell   |           |        0    |
    +--------+-----------+-------------+

    When working with mesh iterators, all entities and connectivity
    are precomputed automatically the first time an iterator is
    created over any given topological dimension or connectivity.

    Note that for efficiency, only entities of dimension zero
    (vertices) and entities of the maximal dimension (cells) exist
    when creating a Mesh. Other entities must be explicitly created
    by calling init(). For example, all edges in a mesh may be created
    by a call to mesh.init(1). Similarly, connectivities such as
    all edges connected to a given vertex must also be explicitly
    created (in this case by a call to mesh.init(0, 1)).

    .. cpp:function:: Mesh(const Mesh& mesh)

        Copy constructor

        *Arguments*
            mesh : a :cpp:class:`Mesh` object.

    .. cpp:function:: Mesh(std::string filename)

        Create mesh from data file.

        *Arguments*
            filename : string, name of file to load. 

    .. cpp:function:: const Mesh& operator=(const Mesh& mesh)

        Assignment

        *Arguments*
            mesh : a :cpp:class:`Mesh` object.

    .. cpp:function:: void all_intersected_entities(const Point & point, uint_set & ids_result) const

        Compute all ids of all cells which are intersected by the given point.

        *Arguments*
            point : a :cpp:class:`Point` object.

            ids_result : set of integers.
                         The cell ids which are intersected are stored in a
                         set for efficiency reasons, to avoid to sort out
                         duplicates later on.

    .. cpp:function:: void all_intersected_entities(const std::vector<Point> & points, uint_set & ids_result) const

        Compute all ids of all cells which are intersected by any point in
        points.

        *Arguments*
            points : vector of :cpp:class:`Point` objects.

            ids_result : set of integers.
                         The cell ids which are intersected are stored in a
                         set for efficiency reasons, to avoid to sort out
                         duplicates later on.

    .. cpp:function:: void all_intersected_entities(const MeshEntity & entity, std::vector<uint> & ids_result) const

        Compute all ids of all cells which are intersected by the given
        entity.

        *Arguments*
            entity : a :cpp:class:`MeshEntity` object.

            ids_result : list of integers.
                         The ids of the intersected cells are saved in a list.
                         This is more efficent than using a set and allows a
                         map between the (external) cell and the intersected
                         cell of the mesh.

    .. cpp:function:: void all_intersected_entities(const std::vector<MeshEntity> & entities, uint_set & ids_result) const

        Compute all id of all cells which are intersected by any entity in the
        vector entities.

        *Arguments*
            entities : vector of :cpp:class:`MeshEntity` objects.

            ids_result : set of integers.
                         The cell ids which are intersected are stored in a
                         set for efficiency reasons, to avoid to sort out
                         duplicates later on.

    .. cpp:function:: void all_intersected_entities(const Mesh & another_mesh, uint_set & ids_result) const

        Compute all ids of all cells which are intersected by another_mesh.

        *Arguments*
            another_mesh : a :cpp:class:`Mesh` object.

            ids_result : set of integers.
                         The cell ids which are intersected are stored in a
                         set for efficiency reasons, to avoid to sort out
                         duplicates later on.

    .. cpp:function:: int any_intersected_entity(const Point & point) const

        Computes only the first id  of the entity, which contains the point.

        *Arguments*
            point : a :cpp:class:`Point` object.

        *Returns*
            integer : the first id of the cell, which contains the point,
            returns -1 if no cell is intersected.

    .. cpp:function:: const uint* cells() const

        Get cell connectivity.

        *Returns*
            array of integers : Connectivity for all cells.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(1,1)
            >>> mesh.coordinates()
            array([[0, 1, 3],
                   [0, 2, 3]])

    .. cpp:function:: void clear()

        Clear all mesh data

    .. cpp:function:: dolfin::uint closest_cell(const Point & point) const

        Computes the index of the cell in the mesh which is closest to the
        point query.

        *Arguments*
            point : a :cpp:class:`Point` object.

        *Returns*
            integer : the index of the cell in the mesh which is closest to
            point.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(1,1)
            >>> point = dolfin.Point(0.0, 2.0)
            >>> mesh.closest_cell(point)
            1

    .. cpp:function:: Point closest_point(const Point & point) const

        Computes the point inside the mesh which is closest to the point
        query.

        *Arguments*
            point : a :cpp:class:`Point` object.

        *Returns*
            :cpp:class:`Point` : the point inside the mesh which is closest to
            the point.

    .. cpp:function:: std::pair<Point,dolfin::uint> closest_point_and_cell(const Point & point) const

        Computes the point inside the mesh and the corresponding cell index
        which are closest to the point query.

        *Arguments*
            point : a :cpp:class:`Point` object.

        *Returns*
            pair <:cpp:class:`Point`, integer> : the point inside the mesh and
            the corresponding cell index which is closest to the point query.

    .. cpp:function:: double* coordinates()

        Get vertex coordinates.

        *Returns*
            array of doubles : Coordinates of all vertices.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(1,1)
            >>> mesh.coordinates()
            array([[ 0.,  0.],
                   [ 1.,  0.],
                   [ 0.,  1.],
                   [ 1.,  1.]])

    .. cpp:function:: const double* coordinates() const

        Return coordinates of all vertices (const version).

    .. cpp:function:: MeshData& data()

        Get mesh data.

        *Returns*
            :cpp:class:`MeshData` : the mesh data object associated with the
            mesh.

    .. cpp:function:: const MeshData& data() const

        Get mesh data (const version).

    .. cpp:function:: MeshGeometry& geometry()

        Get mesh geometry.

        *Returns*
            :cpp:class:`MeshGeometry` : the geometry object associated with the
            mesh.

    .. cpp:function:: const MeshGeometry& geometry() const

        Get mesh geometry (const version).

    .. cpp:function:: double hmax() const

        Compute maximum cell diameter.

        *Returns*
            double : the maximum cell diameter, the diameter is computed as
            two times the circumradius (http://mathworld.wolfram.com).

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.hmax()
            0.70710678118654757

    .. cpp:function:: double hmin() const
        Compute minimum cell diameter.

        *Returns*
            double : the minimum cell diameter, the diameter is computed as
            two times the circumradius (http://mathworld.wolfram.com).

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.hmin()
            0.70710678118654757

    .. cpp:function:: void init() const

        Compute all entities and connectivity.

    .. cpp:function:: uint init(uint dim) const

          Compute entities of given topological dimension.

          *Arguments*
              dim : integer, topological dimension.

          *Returns*
              integer : number of created entities.

    .. cpp:function:: void init(uint d0, uint d1) const

          Compute connectivity between given pair of dimensions.

          *Arguments*
              d0 : integer, topological dimension.

              d1 : integer, topological dimension.

    .. cpp:function:: IntersectionOperator& intersection_operator()

        Get intersectionoperator.

        *Returns*
            :cpp:class:`IntersectionOperator` : the intersection operator
            object associated with the mesh.

    .. cpp:function:: const IntersectionOperator& intersection_operator() const

        Get intersectionoperator (const version).

    .. cpp:function:: void move(BoundaryMesh& boundary, dolfin::ALEType method=hermite)

        Move coordinates of mesh according to new boundary coordinates.

        *Arguments*
            boundary : a :cpp:class:`BoundaryMesh` object.

            method : a :cpp:class:`ALEType` (enum).
                     Method which defines how the coordinates should be moved,
                     default is *hermite*.

    .. cpp:function:: void move(Mesh& mesh, dolfin::ALEType method=hermite)

        Move coordinates of mesh according to adjacent mesh with common global
        vertices.

        *Arguments*
            mesh : a :cpp:class:`Mesh` object.

            method : a :cpp:class:`ALEType` (enum).
                     Method which defines how the coordinates should be moved,
                     default is *hermite*.

    .. cpp:function:: void move(const Function& displacement)

        Move coordinates of mesh according to displacement function. 

        *Arguments*
            function : a :cpp:class:`Function` object.

    .. cpp:function:: uint num_cells() const

        Get number of cells in mesh.

        *Returns*
            integer : number of cells.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_cells()
            8

    .. cpp:function:: uint num_edges() const

        Get number of edges in mesh.

        *Returns*
            integer : number of edges.


        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_edges()
            0
            >>> mesh.init(1)
            16
            >>> mesh.num_edges()
            16

    .. cpp:function:: uint num_entities(uint d) const

        Get number of entities of given topological dimension.

        *Arguments*
            d : integer, topological dimension.

        *Returns*
            integer : number of entities of topological dimension d.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.init(0,1)
            >>> mesh.num_entities(0)
            9
            >>> mesh.num_entities(1)
            16
            >>> mesh.num_entities(2)
            8

    .. cpp:function:: uint num_faces() const

        Get number of faces in mesh.

        *Returns*
            integer : number of faces.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_faces()
            8

    .. cpp:function:: uint num_facets() const

        Get number of facets in mesh.

        *Returns*
            integer : number of facets.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_facets()
            0
            >>> mesh.init(0,1)
            >>> mesh.num_facets()
            16

    .. cpp:function:: uint num_vertices() const

        Get number of vertices in mesh.

        *Returns*
            integer : number of vertices.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.num_vertices()
            9

    .. cpp:function:: void order()

        Order all mesh entities (not needed if "mesh order entities" is set)

        .. seealso::

            UFC documentation (put link here!)

    .. cpp:function:: bool ordered() const

        Check if mesh is ordered.

        *Returns*
            bool : Return true iff topology is ordered according to the UFC
            numbering.

    .. cpp:function:: uint size(uint dim) const

        Get number of entities of given topological dimension.

        *Arguments*
            dim : integer, topological dimension.

        *Returns*
            integer : number of entities of topological dimension d.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.init(0,1)
            >>> mesh.num_entities(0)
            9
            >>> mesh.num_entities(1)
            16
            >>> mesh.num_entities(2)
            8

    .. cpp:function:: void smooth(uint num_iterations=1)

        Smooth internal vertices of mesh by local averaging.

        *Arguments*
            num_iterations : integer,
                             number of iterations to perform smoothing,
                             default value is 1.

    .. cpp:function:: void smooth_boundary(uint num_iterations=1, bool harmonic_smoothing=true)

        Smooth boundary vertices of mesh by local averaging.

        *Arguments*
            num_iterations : integer,
                             number of iterations to perform smoothing,
                             default value is 1.

            harmonic_smoothing : bool,
                                 flag to turn on harmonics smoothing, default
                                 value is true.

    .. cpp:function:: void snap_boundary(const SubDomain& sub_domain, bool harmonic_smoothing=true)

        Snap boundary vertices of mesh to match given sub domain.

        *Arguments*
            sub_domain : a :cpp:class:`SubDomain` object.

            harmonic_smoothing : bool,
                                 a flag to turn on harmonics smoothing,
                                 default value is true.

    .. cpp:function:: std::string str(bool verbose) const

        Informal string representation.

        *Arguments*
            verbose : bool,
                      a flag to turn on additional output.

        *Returns*
            string : an informal representation of the mesh.

        *Example*
            .. warning::

                Not C++ syntax.

            >>> mesh = dolfin.UnitSquare(2,2)
            >>> mesh.str(False)
            '<Mesh of topological dimension 2 (triangles) with 9 vertices and 8 cells, ordered>'

    .. cpp:function:: MeshTopology& topology()

        Get topology associated with mesh.

        *Returns*
            :cpp:class:`MeshTopology` : the topology object associated with the
            mesh.

    .. cpp:function:: const MeshTopology& topology() const

        Get mesh topology (const version).

    .. cpp:function:: CellType& type()

        Get mesh cell type.

        *Returns*
            :cpp:class:`CellType` : the cell type object associated with the
            mesh.

    .. cpp:function:: const CellType& type() const

        Return mesh cell type (const version).

.. .. cpp:function::  ~Mesh() Spinx does not recognize '~' yet!

        Destructor
