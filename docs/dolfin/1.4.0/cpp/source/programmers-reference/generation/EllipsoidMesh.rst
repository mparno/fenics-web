
.. Documentation for the header file dolfin/generation/EllipsoidMesh.h

.. _programmers_reference_cpp_generation_ellipsoidmesh:

EllipsoidMesh.h
===============

.. note::
    
    The documentation on this page was automatically extracted from the
    DOLFIN C++ code and may need to be edited or expanded.
    

.. cpp:class:: EllipsoidMesh

    *Parent class(es)*
    
        * :cpp:class:`Mesh`
        
    Tetrahedral mesh of an ellipsoid. CGAL is used to generate the
    mesh.


    .. cpp:function:: EllipsoidMesh(Point p, std::vector<double> ellipsoid_dims, double cell_size)
    
                EllipsoidMesh mesh(Point(1.0, 2.0, -1.0), dims, 0.2);
        


.. cpp:class:: SphereMesh

    *Parent class(es)*
    
        * :cpp:class:`EllipsoidMesh`
        
    .. cpp:function:: SphereMesh(Point p, double radius, double cell_size)
    
        Create an unstructured :cpp:class:`Mesh` of a sphere
        
        *Arguments*
            center (:cpp:class:`Point`)
                Center of the ellipsoid
            radius (double)
                Axes lengths
            cell_size (double)
                Cell size measure
        
        *Example*
            .. code-block:: c++
        
                // Create sphere with center (1.0, 2.0, -1.0) and
                // radius 3.0. Cell size is 0.2.
                SphereMesh mesh(Point(1.0, 2.0, -1.0), 3.0, 0.2);
        


