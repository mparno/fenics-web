.. Documentation for the header file dolfin/mesh/MeshGeometry.h

.. _programmers_reference_cpp_mesh_meshgeometry:

MeshGeometry.h
==============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and may need to be edited or expanded.

.. cpp:class:: MeshGeometry

    .. cpp:function:: MeshGeometry()
    
        Create empty set of coordinates

    .. cpp:function:: MeshGeometry(const MeshGeometry& geometry)
    
        Copy constructor

    .. cpp:function:: const MeshGeometry& operator= (const MeshGeometry& geometry)
    
        Assignment

    .. cpp:function:: uint dim() const
    
        Return Euclidean dimension of coordinate system

    .. cpp:function:: uint size() const
    
        Return number of coordinates

    .. cpp:function:: double& x(uint n, uint i)
    
        Return value of coordinate n in direction i

    .. cpp:function:: double x(uint n, uint i) const
    
        Return value of coordinate n in direction i

    .. cpp:function:: double* x(uint n)
    
        Return array of values for coordinate n

    .. cpp:function:: const double* x(uint n) const
    
        Return array of values for coordinate n

    .. cpp:function:: double* x()
    
        Return array of values for all coordinates

    .. cpp:function:: const double* x() const
    
        Return array of values for all coordinates

    .. cpp:function:: double* higher_order_x(uint n)
    
        Return array of values for higher order coordinate n

    .. cpp:function:: const double* higher_order_x(uint n) const
    
        Return array of values for higher order coordinate n

    .. cpp:function:: double* higher_order_x()
    
        Return array of values for all higher order coordinates

    .. cpp:function:: const double* higher_order_x() const
    
        Return array of values for all higher order coordinates

    .. cpp:function:: uint num_higher_order_vertices_per_cell() const
    
        Return number of vertices used (per cell) to represent the higher order geometry

    .. cpp:function:: uint* higher_order_cell(uint c)
    
        Return array of higher order vertex indices for a specific higher order cell

    .. cpp:function:: const uint* higher_order_cell(uint c) const
    
        Return array of higher order vertex indices for a specific higher order cell

    .. cpp:function:: uint* higher_order_cells()
    
        Return array of values for all higher order cell data

    .. cpp:function:: const uint* higher_order_cells() const
    
        Return array of values for all higher order cell data

    .. cpp:function:: Point point(uint n) const
    
        Return coordinate n as a 3D point value

    .. cpp:function:: bool* affine_cell_bool()
    
        Return pointer to boolean affine indicator array

    .. cpp:function:: void clear()
    
        Clear all data

    .. cpp:function:: void init(uint dim, uint size)
    
        Initialize coordinate list to given dimension and size

    .. cpp:function:: void init_higher_order_vertices(uint dim, uint size_higher_order)
    
        Initialize higher order coordinate list to given dimension and size

    .. cpp:function:: void init_higher_order_cells(uint num_cells, uint num_dof)
    
        Initialize higher order cell data list to given number of cells and dofs

    .. cpp:function:: void init_affine_indicator(uint num_cells)
    
        Initialize the affine indicator array

    .. cpp:function:: void set_affine_indicator(uint i, bool value)
    
        set affine indicator at index i

    .. cpp:function:: void set(uint n, uint i, double x)
    
        Set value of coordinate n in direction i

    .. cpp:function:: void set_higher_order_coordinates(uint N, uint i, double x)
    
        Set value of higher order coordinate N in direction i

    .. cpp:function:: void set_higher_order_cell_data(uint N, std::vector<uint> vector_cell_data)
    
        Set higher order cell data for cell # N in direction i

    .. cpp:function:: std::string str(bool verbose) const
    
        Return informal string representation (pretty-print)

