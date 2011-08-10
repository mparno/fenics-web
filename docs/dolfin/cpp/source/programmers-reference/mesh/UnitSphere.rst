
.. Documentation for the header file dolfin/mesh/UnitSphere.h

.. _programmers_reference_cpp_mesh_unitsphere:

UnitSphere.h
============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: UnitSphere

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Triangular mesh of the 3D unit sphere.
    
    Given the number of cells (nx, ny, nz) in each direction,
    the total number of tetrahedra will be 6*nx*ny*nz and the
    total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).


