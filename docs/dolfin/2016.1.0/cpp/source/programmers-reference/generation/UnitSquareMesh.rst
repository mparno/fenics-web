
.. Documentation for the header file dolfin/generation/UnitSquareMesh.h

.. _programmers_reference_cpp_generation_unitsquaremesh:

UnitSquareMesh.h
================

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitSquareMesh

    *Parent class(es)*
    
        * :cpp:class:`RectangleMesh`
        
    Triangular mesh of the 2D unit square [0,1] x [0,1].
    Given the number of cells (nx, ny) in each direction,
    the total number of triangles will be 2*nx*ny and the
    total number of vertices will be (nx + 1)*(ny + 1).
    
    std::string diagonal ("left", "right", "right/left", "left/right",
    or "crossed") indicates the direction of the diagonals.


    .. cpp:function:: UnitSquareMesh(std::size_t nx, std::size_t ny, std::string diagonal="right")
    
        Create a uniform finite element :cpp:class:`Mesh` over the unit square
        [0,1] x [0,1].
        
        *Arguments*
            nx (std::size_t)
                Number of cells in horizontal direction.
            ny (std::size_t)
                Number of cells in vertical direction.
            diagonal (std::string)
                Optional argument: A std::string indicating
                the direction of the diagonals.
        
        *Example*
            .. code-block:: c++
        
                UnitSquareMesh mesh1(32, 32);
                UnitSquareMesh mesh2(32, 32, "crossed");
        


    .. cpp:function:: UnitSquareMesh(MPI_Comm comm, std::size_t nx, std::size_t ny, std::string diagonal="right")
    
        Create a uniform finite element :cpp:class:`Mesh` over the unit square
        [0,1] x [0,1].
        
        *Arguments*
            comm (MPI_Comm)
                MPI communicator
            nx (std::size_t)
                Number of cells in horizontal direction.
            ny (std::size_t)
                Number of cells in vertical direction.
            diagonal (std::string)
                Optional argument: A std::string indicating
                the direction of the diagonals.
        
        *Example*
            .. code-block:: c++
        
                UnitSquareMesh mesh1(MPI_COMM_WORLD, 32, 32);
                UnitSquareMesh mesh2(MPI_COMM_WORLD, 32, 32, "crossed");
        


