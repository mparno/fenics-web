
.. Documentation for the header file dolfin/generation/RectangleMesh.h

.. _programmers_reference_cpp_generation_rectanglemesh:

RectangleMesh.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: RectangleMesh

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Triangular mesh of the 2D rectangle spanned by two points p0 and
    p1. Given the number of cells (nx, ny) in each direction, the
    total number of triangles will be 2*nx*ny and the total number
    of vertices will be (nx + 1)*(ny + 1).


    .. cpp:function:: RectangleMesh(const Point& p0, const Point& p1, std::size_t nx, std::size_t ny, std::string diagonal="right")
    
        *Arguments*
            p0 (:cpp:class:`Point`)
                First point.
            p1 (:cpp:class:`Point`)
                Second point.
            nx (double)
                Number of cells in :math:`x`-direction.
            ny (double)
                Number of cells in :math:`y`-direction.
            diagonal (string)
                Direction of diagonals: "left", "right", "left/right", "crossed"
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 8 cells in each direction on the
                // set [-1,2] x [-1,2]
                Point p0(-1, -1);
                Point p1(2, 2);
                RectangleMesh mesh(p0, p1, 8, 8);
        


    .. cpp:function:: RectangleMesh(MPI_Comm comm, const Point& p0, const Point& p1, std::size_t nx, std::size_t ny, std::string diagonal="right")
    
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
            diagonal (string)
                Direction of diagonals: "left", "right", "left/right", "crossed"
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 8 cells in each direction on the
                // set [-1,2] x [-1,2]
                Point p0(-1, -1);
                Point p1(2, 2);
                RectangleMesh mesh(MPI_COMM_WORLD, p0, p1, 8, 8);
        


