
.. Documentation for the header file dolfin/generation/UnitCubeMesh.h

.. _programmers_reference_cpp_generation_unitcubemesh:

UnitCubeMesh.h
==============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitCubeMesh

    *Parent class(es)*
    
        * :cpp:class:`BoxMesh`
        
    Tetrahedral mesh of the 3D unit cube [0,1] x [0,1] x [0,1].
    Given the number of cells (nx, ny, nz) in each direction,
    the total number of tetrahedra will be 6*nx*ny*nz and the
    total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).


    .. cpp:function:: UnitCubeMesh(std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the unit cube
        [0,1] x [0,1] x [0,1].
        
        *Arguments*
            nx (std::size_t)
                Number of cells in :math:`x` direction.
            ny (std::size_t)
                Number of cells in :math:`y` direction.
            nz (std::size_t)
                Number of cells in :math:`z` direction.
        
        *Example*
            .. code-block:: c++
        
                UnitCubeMesh mesh(32, 32, 32);
        


    .. cpp:function:: UnitCubeMesh(MPI_Comm comm, std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the unit cube
        [0,1] x [0,1] x [0,1].
        
        *Arguments*
            comm (MPI_Comm)
                MPI communicator
            nx (std::size_t)
                Number of cells in :math:`x` direction.
            ny (std::size_t)
                Number of cells in :math:`y` direction.
            nz (std::size_t)
                Number of cells in :math:`z` direction.
        
        *Example*
            .. code-block:: c++
        
                UnitCubeMesh mesh(MPI_COMM_WORLD, 32, 32, 32);
        


