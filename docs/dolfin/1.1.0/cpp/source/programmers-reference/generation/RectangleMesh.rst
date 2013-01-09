
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
        
    Triangular mesh of the 2D rectangle (x0, y0) x (x1, y1).
    Given the number of cells (nx, ny) in each direction,
    the total number of triangles will be 2*nx*ny and the
    total number of vertices will be (nx + 1)*(ny + 1).
    
    *Arguments*
        x0 (double)
            :math:`x`-min.
        y0 (double)
            :math:`y`-min.
        x1 (double)
            :math:`x`-max.
        y1 (double)
            :math:`y`-max.
        xn (double)
            Number of cells in :math:`x`-direction.
        yn (double)
            Number of cells in :math:`y`-direction.
        diagonal (string)
            Direction of diagonals: "left", "right", "left/right", "crossed"
    
    *Example*
        .. code-block:: c++
    
            // Mesh with 6 cells in each direction on the
            // set [-1,2] x [-1,2]
            Box mesh(-1, -1, 2, 2, 6, 6;


