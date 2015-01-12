
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
        
    Tetrahedral mesh of the 3D rectangular prism [x0, x1] x [y0, y1]
    x [z0, z1].  Given the number of cells (nx, ny, nz) in each
    direction, the total number of tetrahedra will be 6*nx*ny*nz and
    the total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).


    .. cpp:function:: BoxMesh(double x0, double y0, double z0, double x1, double y1, double z1, std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the rectangular prism
        [x0, x1] x [y0, y1] x [z0, z1].
        
        *Arguments*
            x0 (double)
                :math:`x`-min.
            y0 (double)
                :math:`y`-min.
            z0 (double)
                :math:`z`-min.
            x1 (double)
                :math:`x`-max.
            y1 (double)
                :math:`y`-max.
            z1 (double)
                :math:`z`-max.
            xn (double)
                Number of cells in :math:`x`-direction.
            yn (double)
                Number of cells in :math:`y`-direction.
            zn (double)
                Number of cells in :math:`z`-direction.
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 6 cells in each direction on the
                // set [-1,2] x [-1,2] x [-1,2].
                BoxMesh mesh(-1, -1, -1, 2, 2, 2, 6, 6, 6);
        


    .. cpp:function:: BoxMesh(MPI_Comm comm, double x0, double y0, double z0, double x1, double y1, double z1, std::size_t nx, std::size_t ny, std::size_t nz)
    
        Create a uniform finite element :cpp:class:`Mesh` over the rectangular prism
        [x0, x1] x [y0, y1] x [z0, z1].
        
        *Arguments*
            comm (MPI_Comm)
                MPI communicator
            x0 (double)
                :math:`x`-min.
            y0 (double)
                :math:`y`-min.
            z0 (double)
                :math:`z`-min.
            x1 (double)
                :math:`x`-max.
            y1 (double)
                :math:`y`-max.
            z1 (double)
                :math:`z`-max.
            xn (double)
                Number of cells in :math:`x`-direction.
            yn (double)
                Number of cells in :math:`y`-direction.
            zn (double)
                Number of cells in :math:`z`-direction.
        
        *Example*
            .. code-block:: c++
        
                // Mesh with 6 cells in each direction on the
                // set [-1,2] x [-1,2] x [-1,2].
                BoxMesh mesh(MPI_COMM_WORLD, -1, -1, -1, 2, 2,
                             2, 6, 6, 6);
        


