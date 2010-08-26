.. Documentation for the header file dolfin/mesh/Box.h

.. _programmers_reference_cpp_mesh_box:

Box.h
=====

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: Box

    *Parent class*
    
        * :cpp:class:`Mesh`
        
        Tetrahedral mesh of the 3D  rectangular prism (x0, y0) x (x1, y1) x (x2, y2).
        Given the number of cells (nx, ny, nz) in each direction,
        the total number of tetrahedra will be 6*nx*ny*nz and the
        total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).

