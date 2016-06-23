
.. Documentation for the header file dolfin/mesh/MeshGeometry.h

.. _programmers_reference_cpp_mesh_meshgeometry:

MeshGeometry.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: MeshGeometry

    .. cpp:function:: MeshGeometry()
    
        Create empty set of coordinates


    .. cpp:function:: MeshGeometry(const MeshGeometry& geometry)
    
        Copy constructor


    .. cpp:function:: const MeshGeometry& operator= (const MeshGeometry& geometry)
    
        Assignment


    .. cpp:function:: std::size_t dim() const
    
        Return Euclidean dimension of coordinate system


    .. cpp:function:: std::size_t degree() const
    
        Return polynomial degree of coordinate field


    .. cpp:function:: std::size_t size() const
    
        Return number of coordinates


    .. cpp:function:: std::size_t num_vertices() const
    
        Return the number of vertex coordinates


    .. cpp:function:: std::size_t num_points() const
    
        Return the total number of points in the geometry, located on
        any entity


    .. cpp:function:: double x(std::size_t n, std::size_t i) const
    
        Return value of coordinate with local index n in direction i


    .. cpp:function:: const double* x(std::size_t n) const
    
        Return array of values for coordinate with local index n


    .. cpp:function:: std::vector<double>& x()
    
        Return array of values for all coordinates


    .. cpp:function:: const std::vector<double>& x() const
    
        Return array of values for all coordinates


    .. cpp:function:: Point point(std::size_t n) const
    
        Return coordinate with local index n as a 3D point value


    .. cpp:function:: void init(std::size_t dim, std::size_t degree)
    
        Initialize coordinate list to given dimension and degree


    .. cpp:function:: void init_entities(const std::vector<std::size_t>& num_entities)
    
        Initialise entities. To be called after init


    .. cpp:function:: std::size_t num_entity_coordinates(std::size_t entity_dim) const
    
        Get the number of coordinate points per entity for this degree


    .. cpp:function:: std::size_t get_entity_index(std::size_t entity_dim, std::size_t order, std::size_t index) const
    
        Get the index for an entity point in coordinates


    .. cpp:function:: void set(std::size_t local_index, const double* x)
    
        Set value of coordinate


    .. cpp:function:: std::size_t hash() const
    
        Hash of coordinate values
        
        *Returns*
            std::size_t
                A tree-hashed value of the coordinates over all MPI processes
        


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


