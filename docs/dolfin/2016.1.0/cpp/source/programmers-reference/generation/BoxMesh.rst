
.. Documentation for the header file dolfin/generation/BoxMesh.h

.. _programmers_reference_cpp_generation_boxmesh:

BoxMesh.h
=========

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: BoxMesh

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Tetrahedral mesh of the 3D rectangular prism spanned by two
    points p0 and p1. Given the number of cells (nx, ny, nz) in
    each direction, the total number of tetrahedra will be
    6*nx*ny*nz and the total number of vertices will be
    (nx + 1)*(ny + 1)*(nz + 1).


    .. cpp:function:: BoxMesh(const Point& p0, const Point& p1, std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the rectangular
        prism spanned by the two _Point_s p0 and p1. The order of the
        two points is not important in terms of minimum and maximum
        coordinates.
        
        *Arguments*
            p0 (:cpp:class:`Point`)
                First point.
            p1 (:cpp:class:`Point`)
                Second point.
            nx (double)
                Number of cells in :math:`x`-direction.
            ny (double)
                Number of cells in :math:`y`-direction.
            nz (double)
                Number of cells in :math:`z`-direction.
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 8 cells in each direction on the
                // set [-1,2] x [-1,2] x [-1,2].
                Point p0(-1, -1, -1);
                Point p1(2, 2, 2);
                BoxMesh mesh(p0, p1, 8, 8, 8);
        


    .. cpp:function:: BoxMesh(MPI_Comm comm, const Point& p0, const Point& p1, std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the rectangular
        prism spanned by the two _Point_s p0 and p1. The order of the
        two points is not important in terms of minimum and maximum
        coordinates.
        
        *Arguments*
            comm (MPI_Comm)
                MPI communicator
            p0 (:cpp:class:`Point`)
                First point.
            p1 (:cpp:class:`Point`)
                Second point.
            nx (double)
                Number of cells in :math:`x`-direction.
            ny (double)
                Number of cells in :math:`y`-direction.
            nz (double)
                Number of cells in :math:`z`-direction.
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 8 cells in each direction on the
                // set [-1,2] x [-1,2] x [-1,2].
                Point p0(-1, -1, -1);
                Point p1(2, 2, 2);
                BoxMesh mesh(MPI_COMM_WORLD, p0, p1, 8, 8, 8);
        


