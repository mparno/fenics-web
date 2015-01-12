
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


    .. cpp:function:: std::size_t size() const
    
        Return number of coordinates


    .. cpp:function:: double& x(std::size_t n, std::size_t i)
    
        Return value of coordinate with local index n in direction i


    .. cpp:function:: double x(std::size_t n, std::size_t i) const
    
        Return value of coordinate with local index n in direction i


    .. cpp:function:: double* x(std::size_t n)
    
        Return array of values for coordinate with local index n


    .. cpp:function:: const double* x(std::size_t n) const
    
        Return array of values for coordinate with local index n


    .. cpp:function:: std::vector<double>& x()
    
        Return array of values for all coordinates


    .. cpp:function:: const std::vector<double>& x() const
    
        Return array of values for all coordinates


    .. cpp:function:: Point point(std::size_t n) const
    
        Return coordinate with local index n as a 3D point value


    .. cpp:function:: void clear()
    
        Clear all data


    .. cpp:function:: void init(std::size_t dim, std::size_t size)
    
        Initialize coordinate list to given dimension and size


    .. cpp:function:: std::size_t hash() const
    
        Hash of coordinate values
        
        *Returns*
            std::size_t
                A tree-hashed value of the coordinates over all MPI processes
        


    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)


