.. Documentation for the header file dolfin/mesh/Mesh.h

.. _programmers_reference_cpp_mesh_Mesh:

Mesh.h
======

.. cpp:class:: Mesh

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

    .. cpp:function::  Mesh(const Mesh& mesh)

        Copy constructor

    .. cpp:function::  explicit Mesh(std::string filename)
        
        :param etype: exception type
        :param value: exception value
        :param tb: traceback object
        :param limit: maximum number of stack frames to show
        :type limit: integer or None
        :rtype: list of strings
       
        Create mesh from data file

    .. cpp:function::  const Mesh& operator=(const Mesh& mesh)

        Assignment

    .. cpp:function::  uint num_vertices() const

        Return number of vertices

    .. cpp:function::  uint num_edges() const

        Return number of edges

    .. cpp:function::  uint num_faces() const

        Return number of faces

    .. cpp:function::  uint num_facets() const

        Return number of facets

    .. cpp:function::  uint num_cells() const

        Return number of cells

    .. cpp:function::  uint num_entities(uint d) const

        :arguments:
            * d :type: unsigned int
            * i :type unsigned int:

        :argument d:  unsigned integer (dolfin::uint)
        :param d:  :type: unsigned integer (dolfin::uint)
        :rtype:       unsigned integer
        :returns:     Return number of entities of dimension d

    .. cpp:function::  double* coordinates()

        Return coordinates of all vertices

    .. cpp:function::  const double* coordinates() const

        Return coordinates of all vertices

    .. cpp:function::  const uint* cells() const

        Return connectivity for all cells

    .. cpp:function::  uint size(uint dim) const

        Return number of entities of given topological dimension

    .. cpp:function::  MeshTopology& topology()

        Return mesh topology (non-const version)

    .. cpp:function::  const MeshTopology& topology() const

        Return mesh topology (const version)

    .. cpp:function::  MeshGeometry& geometry()

        Return mesh geometry (non-const version)

    .. cpp:function::  const MeshGeometry& geometry() const

        Return mesh geometry (const version)

    .. cpp:function::  const IntersectionOperator& intersection_operator() const

        Return intersectionoperator (const version)

    .. cpp:function::  IntersectionOperator& intersection_operator()

        Return intersectionoperator (non-const version)

    .. cpp:function::  MeshData& data()

        Return mesh data (non-const version)

    .. cpp:function::  const MeshData& data() const

        Return mesh data (const version)

    .. cpp:function::  CellType& type()

        Return mesh cell type

    .. cpp:function::  const CellType& type() const

        Return mesh cell type

    .. cpp:function::  uint init(uint dim) const

        Compute entities of given topological dimension and return number of entities

    .. cpp:function::  void init(uint d0, uint d1) const

        Compute connectivity between given pair of dimensions

    .. cpp:function::  void init() const

        Compute all entities and connectivity

    .. cpp:function::  void clear()

        Clear all mesh data

    .. cpp:function::  void order()

        Order all mesh entities (not needed if "mesh order entities" is set)

    .. cpp:function::  bool ordered() const

        Return true iff topology is ordered according to the UFC numbering

    .. cpp:function::  void move(BoundaryMesh& boundary, dolfin::ALEType method=hermite)

        Move coordinates of mesh according to new boundary coordinates

    .. cpp:function::  void move(Mesh& mesh, dolfin::ALEType method=hermite)

        Move coordinates of mesh according to adjacent mesh with common global vertices

    .. cpp:function::  void move(const Function& displacement)

        Move coordinates of mesh according to displacement function

    .. cpp:function::  void smooth(uint num_iterations=1)

        Smooth internal vertices of mesh by local averaging

    .. cpp:function::  void smooth_boundary(uint num_iterations=1, bool harmonic_smoothing=true)

        Smooth boundary vertices of mesh by local averaging

    .. cpp:function::  void snap_boundary(const SubDomain& sub_domain, bool harmonic_smoothing=true)

        Snap boundary vertices of mesh to match given sub domain

    .. cpp:function::  void all_intersected_entities(const Point & point, uint_set & ids_result) const

        Compute all id of all cells which are intersects by a \em point.
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.

    .. cpp:function::  void all_intersected_entities(const std::vector<Point> & points, uint_set & ids_result) const

        Compute all id of all cells which are intersects any point in \em points.
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.

    .. cpp:function::  void all_intersected_entities(const MeshEntity & entity, std::vector<uint> & ids_result) const

        Compute all id of all cells which are intersects by a \em entity.
        \param[out] ids_result The ids of the intersected entities are saved in a vector.
        This allows is more efficent than using a set and allows a map between
        the (external) cell and the intersected cell of the mesh. If you
        are only interested in intersection with a list of cells without caring about which
        cell what intersected by which one, use

    .. cpp:function::  void all_intersected_entities(const std::vector<MeshEntity> & entities, uint_set & ids_result) const

        Compute all id of all cells which are intersects by any of the entities in \em entities. This
        \param[out] ids_result The ids of the intersected set are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.

    .. cpp:function::  void all_intersected_entities(const Mesh & another_mesh, uint_set & ids_result) const

        Compute all id of all cells which are intersects by the given mesh \em another_mesh
        \param[out] ids_result The ids of the intersected entities are saved in a set for efficienty
        reasons, to avoid to sort out duplicates later on.

    .. cpp:function::  int any_intersected_entity(const Point & point) const

        Computes only the first id  of the entity, which contains the point. Returns -1 if no cell is intersected.
        @internal @remark This makes the function evaluation significantly faster.

    .. cpp:function::  Point closest_point(const Point & point) const

        Computes the point inside the mesh which are closest to the point query.

    .. cpp:function::  dolfin::uint closest_cell(const Point & point) const

        Computes the index of the cell in the mesh
        which are closest to the point query.

    .. cpp:function::  std::pair<Point,dolfin::uint> closest_point_and_cell(const Point & point) const

        Computes the point inside the mesh and the corresponding cell index
        which are closest to the point query.

    .. cpp:function::  double hmin() const

        Compute minimum cell diameter

    .. cpp:function::  double hmax() const

        Compute maximum cell diameter

    .. cpp:function::  std::string str(bool verbose) const

        Return informal string representation (pretty-print)

.. .. cpp:function::  ~Mesh() Spinx does not recognize '~' yet!

        Destructor
