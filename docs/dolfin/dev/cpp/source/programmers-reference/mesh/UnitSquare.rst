
.. Documentation for the header file dolfin/mesh/UnitSquare.h

.. _programmers_reference_cpp_mesh_unitsquare:

UnitSquare.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitSquare

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Triangular mesh of the 2D unit square [0,1] x [0,1].
    Given the number of cells (nx, ny) in each direction,
    the total number of triangles will be 2*nx*ny and the
    total number of vertices will be (nx + 1)*(ny + 1).
    
    std::string diagonal ("left", "right", "right/left", "left/right",
    or "crossed") indicates the direction of the diagonals.


    .. cpp:function:: UnitSquare(uint nx, uint ny, std::string diagonal="right")
    
        Define a uniform finite element :cpp:class:`Mesh` over the unit square
        [0,1] x [0,1].
        
        *Arguments*
            nx (uint)
                Number of cells in horizontal direction.
            ny (uint)
                Number of cells in vertical direction.
            diagonal (std::string)
                Optional argument: A std::string indicating
                the direction of the diagonals.
        
        *Example*
            .. code-block:: c++
        
                UnitSquare mesh1(32, 32);
                UnitSquare mesh2(32, 32, "crossed");
        


