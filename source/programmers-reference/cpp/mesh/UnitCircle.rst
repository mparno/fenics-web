.. Documentation for the header file dolfin/mesh/UnitCircle.h

.. _programmers_reference_cpp_mesh_Mesh:

UnitCircle.h
============

.. note::

    The documentation on this page was automatically extracted from
    the DOLFIN C++ code and needs to be edited and expanded.

.. cpp:class:: UnitCircle

    *Parent class*
    
        * :cpp:class:`Mesh`
        
        Triangular mesh of the 2D unit circle.
        Given the number of cells (nx, ny) in each direction,
        the total number of triangles will be 2*nx*ny and the
        total number of vertices will be (nx + 1)*(ny + 1).
        std::string diagonal ("left", "right" or "crossed") indicates the
        direction of the diagonals.
        std:string transformation ("maxn", "sumn" or "rotsumn")

